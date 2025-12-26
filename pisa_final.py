#!/usr/bin/env python3
# pisa.py
#
# Batch-runner for CCP4 PISA interface analysis on all *.pdb files
# in the current directory.
#
# Configuration:
#   - No settings file support; all configuration via CLI.
#
# Added outputs (contacts.csv):
#   - Interface energetics (best interface involving requested chain): ΔG of dissociation, solvation energy gain (if present in PISA XML)
#   - Interface area in Å² (sum across interfaces involving requested chain, if present in PISA XML; otherwise blank)
#   - Interface residue count and % polar/hydrophobic (based on chain residues with BSA>0)
#
# Outputs created:
#   - contacts.csv
#   - pisa_analysis.log
#   - pisa_xml_files/  (directory containing all generated XML files)

import os
import sys
import shutil
import subprocess
import argparse
import math
import logging
import xml.etree.ElementTree as ET
from concurrent.futures import ThreadPoolExecutor, as_completed

HYDROPHOBIC = {"ALA","VAL","LEU","ILE","MET","PHE","TRP","PRO","GLY"}
POLAR_UNCHARGED = {"SER","THR","ASN","GLN","TYR","CYS"}
CHARGED = {"ASP","GLU","LYS","ARG","HIS"}

def setup_logging():
    log_file = "pisa_analysis.log"
    logging.basicConfig(
        level=logging.INFO,
        format="%(asctime)s [%(levelname)s] %(message)s",
        handlers=[logging.FileHandler(log_file), logging.StreamHandler(sys.stdout)],
    )
    logging.info("Logging is set up.")

def parse_arguments():
    parser = argparse.ArgumentParser(description="Run PISA on all .pdb files in the current directory.")
    parser.add_argument("--chain", "-C", default="A", help="Chain ID to analyze (default: A). Example: --chain H")
    parser.add_argument("--residue_counter", "-r", type=int, default=0,
                        help="Only count contacts for residues with seq_num >= this value (default: 0)")
    parser.add_argument("--pisa", default=None, help="Explicit path to PISA executable (highest priority).")
    parser.add_argument("--use_module", action="store_true", help="Attempt module load (requires --ccp4_module).")
    parser.add_argument("--ccp4_module", default=None,
                        help="CCP4 module name to load (e.g. ccp4/8.0.010). Only used with --use_module.")
    parser.add_argument("--max_workers", type=int, default=None, help="Override max_workers.")
    parser.add_argument("--preflight_only", action="store_true", help="Resolve PISA executable and exit.")
    return parser.parse_args()

def _get_physical_cores_fallback():
    try:
        import psutil  # optional dependency
        n = psutil.cpu_count(logical=False)
        if n:
            return int(n)
    except Exception:
        pass
    return max(1, (os.cpu_count() or 1))

def _get_available_mem_gb_linux_fallback():
    try:
        with open("/proc/meminfo", "r") as f:
            kv = {}
            for line in f:
                if ":" not in line:
                    continue
                k, rest = line.split(":", 1)
                parts = rest.strip().split()
                if not parts:
                    continue
                kv[k.strip()] = int(parts[0])  # kB
        if "MemAvailable" in kv:
            return kv["MemAvailable"] / (1024 * 1024)
        if "MemTotal" in kv:
            return kv["MemTotal"] / (1024 * 1024)
    except Exception:
        pass
    return 4.0

def estimate_max_workers():
    physical = _get_physical_cores_fallback()
    avail_gb = _get_available_mem_gb_linux_fallback()
    mem_limited = max(1, int(avail_gb // 1.5))
    est = min(physical, mem_limited, 32)
    return max(1, est)

def load_environment(*, ccp4_module, pisa_path, use_module):
    if pisa_path:
        pisa_path = os.path.abspath(pisa_path)
        if not (os.path.isfile(pisa_path) and os.access(pisa_path, os.X_OK)):
            raise RuntimeError(f"PISA path provided but is not an executable file: {pisa_path}")
        os.environ["PISA_EXE"] = pisa_path
        logging.info(f"Using explicit PISA executable: {pisa_path}")
        return

    if use_module:
        if not ccp4_module:
            raise RuntimeError("--use_module was set but --ccp4_module was not provided.")
        logging.info(f"Attempting to load CCP4 module: {ccp4_module}")
        try:
            result = subprocess.run(
                f"module load {ccp4_module} >/dev/null 2>&1 && env",
                shell=True,
                check=True,
                stdout=subprocess.PIPE,
                executable="/bin/bash",
            )
            env_vars = {}
            for line in result.stdout.decode("utf-8").splitlines():
                if "=" in line:
                    k, v = line.split("=", 1)
                    env_vars[k] = v
            os.environ.update(env_vars)
            logging.info("Module environment loaded into current process.")
        except subprocess.CalledProcessError as e:
            raise RuntimeError(
                f"Failed to load module '{ccp4_module}'. "
                f"Either fix the module name, set --pisa /path/to/pisa, or ensure pisa is already on PATH."
            ) from e
    else:
        logging.info("Not loading any module (default). Relying on PATH for PISA.")

def preflight_check_pisa():
    pisa_exe = os.environ.get("PISA_EXE")
    if pisa_exe:
        if not (os.path.isfile(pisa_exe) and os.access(pisa_exe, os.X_OK)):
            raise RuntimeError(f"PISA_EXE is set but not executable: {pisa_exe}")
        return pisa_exe

    found = shutil.which("pisa")
    if not found:
        raise RuntimeError(
            "PISA executable not found. Remediation options:\n"
            "  - Put 'pisa' on PATH (module load, etc.)\n"
            "  - OR run with --pisa /full/path/to/pisa\n"
            "  - OR run with --use_module --ccp4_module <name>\n"
        )
    return found

def _assert_chain_present(xml_file, requested_chain):
    try:
        root = ET.parse(xml_file).getroot()
    except Exception as e:
        raise RuntimeError(f"Cannot parse XML for chain validation: {xml_file}. Error: {e}") from e

    chains = set()
    for molecule in root.findall(".//molecule"):
        cid = molecule.findtext("chain_id", default="").strip()
        if cid:
            chains.add(cid)

    if requested_chain not in chains:
        chains_sorted = ", ".join(sorted(chains)) if chains else "(none found)"
        raise RuntimeError(
            f"Requested chain '{requested_chain}' not found in {xml_file}. Chains present: {chains_sorted}."
        )

def parse_xml_residues_and_bonds(xml_file, chain_id, residue_counter):
    try:
        root = ET.parse(xml_file).getroot()
    except Exception as e:
        logging.error(f"Error parsing {xml_file}: {e}")
        return [], 0, 0

    residues = []
    h_bonds_count = 0
    salt_bridges_count = 0

    for molecule in root.findall(".//molecule"):
        current_chain_id = (molecule.findtext("chain_id", default="") or "").strip()
        for residue in molecule.findall(".//residue"):
            res_name = (residue.findtext("name", default="") or "").strip()
            seq_num = int(residue.findtext("seq_num", default="0") or "0")
            asa = float(residue.findtext("asa", default="0") or "0")
            bsa = float(residue.findtext("bsa", default="0") or "0")
            solv_en = float(residue.findtext("solv_en", default="0") or "0")
            residues.append((current_chain_id, res_name, seq_num, asa, bsa, solv_en))

    for interface in root.findall(".//interface"):
        hb = interface.find("h-bonds")
        if hb is not None:
            for bond in hb.findall("bond"):
                c1 = (bond.findtext("chain-1", default="") or "").strip()
                c2 = (bond.findtext("chain-2", default="") or "").strip()
                if c1 == chain_id or c2 == chain_id:
                    s1 = int(bond.findtext("seqnum-1", default="0") or "0")
                    s2 = int(bond.findtext("seqnum-2", default="0") or "0")
                    if c1 == chain_id and s1 >= residue_counter:
                        h_bonds_count += 1
                    if c2 == chain_id and s2 >= residue_counter:
                        h_bonds_count += 1

        sb = interface.find("salt-bridges")
        if sb is not None:
            for bond in sb.findall("bond"):
                c1 = (bond.findtext("chain-1", default="") or "").strip()
                c2 = (bond.findtext("chain-2", default="") or "").strip()
                if c1 == chain_id or c2 == chain_id:
                    s1 = int(bond.findtext("seqnum-1", default="0") or "0")
                    s2 = int(bond.findtext("seqnum-2", default="0") or "0")
                    if c1 == chain_id and s1 >= residue_counter:
                        salt_bridges_count += 1
                    if c2 == chain_id and s2 >= residue_counter:
                        salt_bridges_count += 1

    return residues, h_bonds_count, salt_bridges_count

def calculate_total_bsa_score(residues, chain_id, residue_counter):
    total_bars = 0
    for current_chain_id, _, seq_num, asa, bsa, _ in residues:
        if current_chain_id == chain_id and seq_num >= residue_counter and asa != 0:
            pct = (bsa / asa) * 100.0
            total_bars += int(pct / 10.0)
    return total_bars

def compute_interface_residue_stats(residues, chain_id, residue_counter):
    iface = []
    for current_chain_id, res_name, seq_num, asa, bsa, solv_en in residues:
        if current_chain_id != chain_id:
            continue
        if seq_num < residue_counter:
            continue
        if bsa and bsa > 0:
            iface.append(res_name.strip().upper())

    n = len(iface)
    if n == 0:
        return 0, "", ""

    hyd = sum(1 for r in iface if r in HYDROPHOBIC)
    polar = sum(1 for r in iface if (r in POLAR_UNCHARGED or r in CHARGED))
    pct_hyd = 100.0 * hyd / n
    pct_polar = 100.0 * polar / n
    return n, f"{pct_polar:.2f}", f"{pct_hyd:.2f}"

def _float_or_none(x):
    try:
        if x is None:
            return None
        s = str(x).strip()
        if s == "":
            return None
        return float(s)
    except Exception:
        return None

def _find_first_float_in_interface(interface_elem, candidate_tags):
    for tag in candidate_tags:
        val = interface_elem.findtext(tag)
        f = _float_or_none(val)
        if f is not None:
            return f
    for child in list(interface_elem):
        t = (child.tag or "").lower()
        for tag in candidate_tags:
            if t == tag.lower():
                f = _float_or_none(child.text)
                if f is not None:
                    return f
    return None

def parse_interface_energetics_and_area(xml_file, chain_id):
    """
    Extract interface area and energetics for interfaces involving the requested chain.

    XML schema (your PISA build):
      - interface area: <int_area>
      - solvation term: <int_solv_en>
      - stabilisation energy: <stab_en>
      - overlap: <overlap>
      - pvalue: <pvalue>
      - participating chains: <interface><molecule><chain_id>...</chain_id></molecule>...

    We select the "best" interface as the most stabilising (most negative stab_en).

    Returns (6-tuple):
      iface_count,
      total_area_A2_str,
      dg_dissociation_str    (RAW <stab_en> for best interface),
      solvation_gain_str     (RAW <int_solv_en> for best interface),
      overlap_str            (RAW <overlap> for best interface),
      specificity_str        (-log10(pvalue) for best interface; blank if missing/invalid)
    """
    try:
        root = ET.parse(xml_file).getroot()
    except Exception as e:
        logging.error(f"Error parsing {xml_file} for interface energetics: {e}")
        return 0, "", "", "", "", ""

    iface_count = 0
    total_area = 0.0
    have_area = False

    best_stab_en = None
    best_solv_en = None
    best_overlap = None
    best_pvalue = None

    for interface in root.findall(".//interface"):
        chains = set()
        for m in interface.findall("molecule"):
            cid = (m.findtext("chain_id") or "").strip()
            if cid:
                chains.add(cid)

        if chain_id not in chains:
            continue

        iface_count += 1

        int_area = _float_or_none(interface.findtext("int_area"))
        if int_area is not None:
            total_area += int_area
            have_area = True

        stab_en = _float_or_none(interface.findtext("stab_en"))
        solv_en = _float_or_none(interface.findtext("int_solv_en"))
        overlap_txt = (interface.findtext("overlap") or "").strip()
        # PISA reports overlap as Yes/No (string) in many builds
        overlap = overlap_txt if overlap_txt else ""
        pvalue = _float_or_none(interface.findtext("pvalue"))

        if stab_en is not None and (best_stab_en is None or stab_en < best_stab_en):
            best_stab_en = stab_en
            best_solv_en = solv_en
            best_overlap = overlap
            best_pvalue = pvalue

    total_area_str = f"{total_area:.2f}" if have_area else ""
    dg_str = f"{best_stab_en:.2f}" if best_stab_en is not None else ""
    solv_str = f"{best_solv_en:.3f}" if best_solv_en is not None else ""
    overlap_str = str(best_overlap) if best_overlap not in (None, "") else ""

    specificity_str = ""
    if best_pvalue is not None and best_pvalue > 0:
        try:
            spec = -math.log10(best_pvalue)
            if math.isfinite(spec):
                specificity_str = f"{spec:.3f}"
        except Exception:
            specificity_str = ""

    return iface_count, total_area_str, dg_str, solv_str, overlap_str, specificity_str

def process_pdb_file(pdb_file, chain_id, residue_counter):
    base_filename = os.path.splitext(pdb_file)[0]
    pisa_exe = os.environ.get("PISA_EXE") or "pisa"
    try:
        logging.info(f"Starting PISA analysis for {pdb_file}")
        subprocess.run([pisa_exe, pdb_file, "-analyse", pdb_file],
                       stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL, check=True)

        xml_output = f"{base_filename}.xml"
        with open(xml_output, "w") as xml_f:
            subprocess.run([pisa_exe, pdb_file, "-xml", "interfaces"],
                           stdout=xml_f, stderr=subprocess.DEVNULL, check=True)

        _assert_chain_present(xml_output, chain_id)

        subprocess.run([pisa_exe, pdb_file, "-erase"],
                       stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL, check=True)
    except Exception as e:
        logging.error(f"Failed processing {pdb_file}: {e}")
        return base_filename, 0, 0, 0, 0, "", "", "", "", "", 0, "", ""

    residues, h_bonds_count, salt_bridges_count = parse_xml_residues_and_bonds(xml_output, chain_id, residue_counter)
    bsa_score = calculate_total_bsa_score(residues, chain_id, residue_counter)
    iface_res_n, pct_polar, pct_hydrophobic = compute_interface_residue_stats(residues, chain_id, residue_counter)
    iface_count, iface_area_A2, dg_diss, solv_gain, overlap, specificity = parse_interface_energetics_and_area(xml_output, chain_id)

    return (base_filename, bsa_score, salt_bridges_count, h_bonds_count,
            iface_count, iface_area_A2, dg_diss, solv_gain, overlap, specificity,
            iface_res_n, pct_polar, pct_hydrophobic)

def main():
    setup_logging()
    args = parse_arguments()

    load_environment(ccp4_module=args.ccp4_module, pisa_path=args.pisa, use_module=bool(args.use_module))
    pisa_resolved = preflight_check_pisa()
    logging.info(f"Using PISA executable: {pisa_resolved}")

    if args.preflight_only:
        logging.info("--preflight_only requested; exiting.")
        return

    max_workers = int(args.max_workers) if args.max_workers is not None else estimate_max_workers()
    max_workers = max(1, max_workers)
    logging.info(f"max_workers={max_workers}")

    chain_id = (args.chain or "A").strip()
    residue_counter = max(0, int(args.residue_counter or 0))

    pdb_files = [f for f in os.listdir(".") if f.endswith(".pdb")]
    total_files = len(pdb_files)
    logging.info(f"Found {total_files} .pdb files in '{os.getcwd()}'.")

    if total_files == 0:
        logging.warning("No .pdb files found. Exiting.")
        return

    results = []
    processed = 0
    last_decile = -1

    with ThreadPoolExecutor(max_workers=max_workers) as executor:
        futures = {executor.submit(process_pdb_file, f, chain_id, residue_counter): f for f in pdb_files}
        for future in as_completed(futures):
            res = future.result()
            results.append(res)
            processed += 1
            pct = (processed / total_files) * 100.0
            decile = int(pct // 10)
            if decile != last_decile:
                last_decile = decile
                logging.info(f"Processed {processed}/{total_files} files ({int(pct)}% complete)")

    xml_dir = "pisa_xml_files"
    os.makedirs(xml_dir, exist_ok=True)
    for f in os.listdir("."):
        if f.endswith(".xml"):
            try:
                shutil.move(f, os.path.join(xml_dir, f))
            except Exception as e:
                logging.error(f"Failed to move {f}: {e}")
    with open("contacts.csv", "w") as out:
        out.write(
            "binder,bsa_score,salt_bridges,h_bonds,"
            "interface_count,interface_area_A2,"
            "dg_dissociation,solvation_energy_gain,"
            "overlap,specificity,"
            "interface_residue_count,pct_polar,pct_hydrophobic\n"
        )
        for (base_name, bsa_score, salt_bridges, h_bonds,
             iface_count, iface_area_A2, dg_diss, solv_gain, overlap, specificity,
             iface_res_n, pct_polar, pct_hydrophobic) in results:
            out.write(
                f"{base_name},{bsa_score},{salt_bridges},{h_bonds},"
                f"{iface_count},{iface_area_A2},"
                f"{dg_diss},{solv_gain},"
                f"{overlap},{specificity},"
                f"{iface_res_n},{pct_polar},{pct_hydrophobic}\n"
            )
    logging.info("Wrote contacts.csv")

if __name__ == "__main__":
    main()


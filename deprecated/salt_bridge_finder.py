import os
import xml.etree.ElementTree as ET
import logging
import csv
import argparse
import re
from concurrent.futures import ProcessPoolExecutor

# Configure logging
logging.basicConfig(
    filename='salt_bridge_filter.log',
    filemode='w',
    level=logging.INFO,
    format='%(levelname)s:%(message)s'
)

def process_single_file(filepath, target_chain_A, target_chain_B):
    """
    Processes a single XML file and extracts salt bridges.

    Parameters:
        filepath (str): Path to the XML file.
        target_chain_A (str): Chain ID for Chain A.
        target_chain_B (str): Chain ID for Chain B.

    Returns:
        list: List of salt bridge data in the format:
              [formatted_pdb_id, Residue1, Residue2, Distance]
    """
    results = []
    filename = os.path.basename(filepath)
    pdb_id = filename.replace('.xml', '').lower()

    try:
        tree = ET.parse(filepath)
        root = tree.getroot()
    except ET.ParseError as e:
        logging.error(f"Failed to parse {filename}: {e}. Skipping file.")
        return results

    interfaces = root.findall('.//interface')
    if not interfaces:
        logging.warning(f"No interfaces found in {filename}. Skipping.")
        return results

    for interface in interfaces:
        molecules = interface.findall('molecule')
        chain_residues = {}

        for molecule in molecules:
            chain_id_elem = molecule.find('chain_id')
            if chain_id_elem is None or chain_id_elem.text is None:
                logging.warning(f"Missing <chain_id> in {filename}. Skipping molecule.")
                continue

            chain_id = chain_id_elem.text.strip()
            residues = molecule.findall('.//residues/residue')
            residue_seq_nums = []

            for res in residues:
                seq_num_elem = res.find('seq_num')
                if seq_num_elem is None or seq_num_elem.text is None:
                    continue
                try:
                    seq_num = int(seq_num_elem.text.strip())
                    residue_seq_nums.append(seq_num)
                except ValueError:
                    continue

            if residue_seq_nums:
                chain_residues[chain_id] = residue_seq_nums

        if target_chain_B not in chain_residues:
            logging.warning(f"Chain {target_chain_B} not found in {filename}. Skipping.")
            continue

        first_seq_num_B = min(chain_residues[target_chain_B])
        target_residue_B = first_seq_num_B + 30

        salt_bridges = interface.find('salt-bridges')
        if salt_bridges is None:
            continue

        for bond in salt_bridges.findall('bond'):
            chain1_elem = bond.find('chain-1')
            res1_elem = bond.find('res-1')
            seqnum1_elem = bond.find('seqnum-1')

            chain2_elem = bond.find('chain-2')
            res2_elem = bond.find('res-2')
            seqnum2_elem = bond.find('seqnum-2')

            dist_elem = bond.find('dist')

            if None in (chain1_elem, res1_elem, seqnum1_elem, chain2_elem, res2_elem, seqnum2_elem, dist_elem):
                continue

            try:
                seqnum1 = int(seqnum1_elem.text.strip())
                seqnum2 = int(seqnum2_elem.text.strip())
                dist = float(dist_elem.text.strip())
            except ValueError:
                continue

            if (chain1_elem.text.strip() == target_chain_B and seqnum1 == target_residue_B and chain2_elem.text.strip() == target_chain_A) or \
               (chain2_elem.text.strip() == target_chain_B and seqnum2 == target_residue_B and chain1_elem.text.strip() == target_chain_A):
                res1 = f"{chain1_elem.text.strip()}/{res1_elem.text.strip()}{seqnum1}"
                res2 = f"{chain2_elem.text.strip()}/{res2_elem.text.strip()}{seqnum2}"
                dist_rounded = round(dist, 2)

                formatted_pdb_id = pdb_id # Already lowercase
                logging.info(f"{formatted_pdb_id}: {res1} ↔ {res2} at {dist_rounded} Å")

                results.append([formatted_pdb_id, res1, res2, dist_rounded])

    return results

def read_compute_settings(settings_file):
    """
    Reads compute settings from a file.

    Parameters:
        settings_file (str): Path to the settings file.

    Returns:
        dict: A dictionary containing compute settings (e.g., max_workers, batch_size).
    """
    settings = {'max_workers': 4, 'batch_size': None}  # Defaults
    if os.path.exists(settings_file):
        with open(settings_file, 'r') as f:
            lines = f.readlines()
        # Only parse lines under "Recommended settings:"
        start_parsing = False
        for line in lines:
            if "Recommended settings:" in line:
                start_parsing = True
                continue
            if start_parsing:
                if 'max_workers' in line:
                    settings['max_workers'] = int(line.split(':')[1].strip())
                if 'batch_size' in line:
                    settings['batch_size'] = int(line.split(':')[1].strip())
    return settings

def find_pdb_files_with_specific_salt_bridge(xml_directory, target_chain_A, target_chain_B, output_csv, max_workers):
    """
    Processes XML files in parallel to extract salt bridges.

    Parameters:
        xml_directory (str): Path to the directory containing XML files.
        target_chain_A (str): Chain ID for Chain A.
        target_chain_B (str): Chain ID for Chain B.
        output_csv (str): Path to the output CSV file.
        max_workers (int): Number of parallel workers (processes).
    """
    results = []

    # List all XML files in the directory
    filepaths = [os.path.join(xml_directory, f) for f in os.listdir(xml_directory) if f.endswith('.xml')]

    with ProcessPoolExecutor(max_workers=max_workers) as executor:
        futures = [executor.submit(process_single_file, filepath, target_chain_A, target_chain_B) for filepath in filepaths]

        for future in futures:
            try:
                results.extend(future.result())
            except Exception as e:
                logging.error(f"Error during processing: {e}")

    # Write results to CSV
    if results:
        with open(output_csv, mode='w', newline='') as csvfile:
            csv_writer = csv.writer(csvfile)
            csv_writer.writerow(['PDB_ID', 'Residue1', 'Residue2', 'Distance (Å)'])
            csv_writer.writerows(results)
        print(f"Results written to {output_csv}")
    else:
        print("No qualifying salt bridges found.")

def main():
    parser = argparse.ArgumentParser(description="Filter PISA XML files for specific salt bridges.")
    parser.add_argument('xml_directory', type=str, nargs='?', help="Path to the directory containing PISA XML files.")
    parser.add_argument('--output_csv', type=str, default='qualifying_salt_bridges.csv',
                        help="Path to the output CSV file (default: qualifying_salt_bridges.csv).")
    parser.add_argument('--chain_a', type=str, default='A', help="Chain ID for Chain A (default: A).")
    parser.add_argument('--chain_b', type=str, default='B', help="Chain ID for Chain B (default: B).")
    parser.add_argument('--settings_file', type=str, default='compute_settings.txt',
                        help="Path to the compute settings file (default: compute_settings.txt).")

    args = parser.parse_args()

    if not args.xml_directory:
        parser.print_help()
        return

    if not os.path.isdir(args.xml_directory):
        logging.critical(f"The directory {args.xml_directory} does not exist.")
        print(f"Error: The directory {args.xml_directory} does not exist. Check the path and try again.")
        return

    compute_settings = read_compute_settings(args.settings_file)
    max_workers = compute_settings['max_workers']

    find_pdb_files_with_specific_salt_bridge(
        args.xml_directory,
        args.chain_a,
        args.chain_b,
        args.output_csv,
        max_workers
    )

if __name__ == "__main__":
    main()


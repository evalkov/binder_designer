#!/usr/bin/env python3

import os
import shutil
import subprocess
import sys
import csv
import argparse
from pathlib import Path
from concurrent.futures import ThreadPoolExecutor, as_completed
import logging
from typing import Optional, Tuple, List, Dict

# Constants
AF2_PAE_INTERACT = 10

# Configure logging
logging.basicConfig(
    level=logging.INFO,
    format='%(asctime)s [%(levelname)s] %(message)s',
    handlers=[
        logging.FileHandler("process_binders.log"),
        logging.StreamHandler()
    ]
)

# Mapping from three-letter to one-letter amino acid codes
THREE_TO_ONE = {
    "ALA": "A", "CYS": "C", "ASP": "D", "GLU": "E", "PHE": "F", "GLY": "G",
    "HIS": "H", "ILE": "I", "LYS": "K", "LEU": "L", "MET": "M", "ASN": "N",
    "PRO": "P", "GLN": "Q", "ARG": "R", "SER": "S", "THR": "T", "VAL": "V",
    "TRP": "W", "TYR": "Y",
}

def load_python_module():
    """
    Load the Python module and determine the Python executable.
    Exits the script if loading fails.
    """
    command = "module load python"
    logging.info(f"Loading Python module with: {command}")
    result = subprocess.run(command, shell=True)
    if result.returncode != 0:
        logging.error("Failed to load Python module. Exiting.")
        sys.exit(1)

    # Verify the Python executable is available
    python_exec = subprocess.run(
        "which python", shell=True, stdout=subprocess.PIPE, universal_newlines=True
    ).stdout.strip()
    if not python_exec:
        logging.warning("Python executable not found after loading module. Trying python3.")
        python_exec = subprocess.run(
            "which python3", shell=True, stdout=subprocess.PIPE, universal_newlines=True
        ).stdout.strip()

    if not python_exec:
        logging.error("No valid Python executable found. Exiting.")
        sys.exit(1)

    logging.info(f"Using Python executable: {python_exec}")
    return python_exec

def execute_command(command: str, cwd: Optional[str] = None):
    """
    Execute a shell command, optionally changing the working directory.
    Exits the script if the command fails.

    Args:
        command (str): The shell command to execute.
        cwd (Optional[str]): The working directory to execute the command in.
    """
    if cwd:
        logging.info(f"Changing directory to {cwd}")
    logging.info(f"Executing command: {command}")
    result = subprocess.run(command, shell=True, cwd=cwd)
    if result.returncode != 0:
        logging.error(f"Command failed: {command}")
        sys.exit(1)

def combine_sc_files(sc_directory: Path) -> List[str]:
    """
    Combine all out_*.sc files into a single list of lines, with special handling for out_1.sc.

    Args:
        sc_directory (Path): Directory containing out_*.sc files.

    Returns:
        List[str]: Combined and filtered lines from all out_*.sc files.
    """
    out_sc = sc_directory / "out.sc"
    out_1_sc = sc_directory / "out_1.sc"

    combined_lines = []
    if out_1_sc.exists():
        with out_1_sc.open("r") as infile:
            header = infile.readline().strip()
            combined_lines.append(header)
            for line in infile:
                if line.startswith("SCORE:"):
                    combined_lines.append(line)
        # Process other out_*.sc files
        for file in sorted(f for f in sc_directory.iterdir() if f.name.startswith("out_") and f.name != "out_1.sc" and f.suffix == ".sc"):
            with file.open("r") as infile:
                infile.readline()  # Skip header
                for line in infile:
                    if line.startswith("SCORE:"):
                        combined_lines.append(line)
        logging.info(f"Combined {len(combined_lines)} lines from out_*.sc files.")
    else:
        logging.error("out_1.sc doesn't exist. Exiting.")
        sys.exit(1)

    return combined_lines

def filter_sc_lines(combined_lines: List[str]) -> List[List[str]]:
    """
    Filter combined .sc lines based on pae_interaction and extract required fields.

    Args:
        combined_lines (List[str]): Combined lines from out_*.sc files.

    Returns:
        List[List[str]]: Filtered data with required fields.
    """
    if not combined_lines:
        logging.error("No lines to process in combined_sc_files.")
        sys.exit(1)

    # Extract and validate header
    header_line = combined_lines[0]
    header = header_line.replace("SCORE:", "").strip().split()
    target_columns = [
        "binder_aligned_rmsd",
        "pae_binder",
        "pae_interaction",
        "plddt_binder",
        "description",
    ]
    missing_columns = [col for col in target_columns if col not in header]

    if missing_columns:
        logging.error(f"Missing expected columns in the header: {missing_columns}")
        sys.exit(1)

    indices = [header.index(col) for col in target_columns]
    filtered_data = []

    for line in combined_lines[1:]:
        if not line.startswith("SCORE:"):
            continue
        columns = line.replace("SCORE:", "").strip().split()
        try:
            pae_interaction = float(columns[indices[2]])  # "pae_interaction"
            if pae_interaction < AF2_PAE_INTERACT:
                selected_data = [columns[i] for i in indices]
                filtered_data.append(selected_data)
        except (ValueError, IndexError) as e:
            logging.warning(f"Skipping line due to error: {e}")
            continue

    logging.info(f"Filtered down to {len(filtered_data)} entries with pae_interaction < {AF2_PAE_INTERACT}.")
    return filtered_data

def parse_settings_file(settings_file: Path) -> Dict[str, str]:
    """
    Parse the compute_settings.txt file to extract configurations.

    Args:
        settings_file (Path): Path to the settings file.

    Returns:
        Dict[str, str]: A dictionary containing configuration parameters.
    """
    if not settings_file.is_file():
        logging.error(f"Settings file does not exist: {settings_file}")
        raise FileNotFoundError(f"Settings file does not exist: {settings_file}")

    settings = {}
    with settings_file.open('r') as f:
        for line in f:
            line = line.strip()
            if not line or line.startswith('#'):
                continue  # Skip empty lines and comments
            if ':' in line:
                key, value = line.split(':', 1)
                settings[key.strip()] = value.strip()

    return settings

def extract_sequence(filename: Path, chain: str = "A") -> Optional[str]:
    """
    Extract the amino acid sequence for a specific chain from a PDB file.

    Args:
        filename (Path): Path to the PDB file.
        chain (str): Chain identifier to extract the sequence from.

    Returns:
        Optional[str]: The extracted amino acid sequence or None if not found.
    """
    sequence = []
    try:
        with filename.open("r") as pdb_file:
            for line in pdb_file:
                if line.startswith("ATOM"):
                    atom_name = line[12:16].strip()
                    chain_id = line[21].strip()
                    if atom_name == "CA" and chain_id == chain:
                        residue = line[17:20].strip()
                        amino_acid = THREE_TO_ONE.get(residue, "X")
                        sequence.append(amino_acid)
        if sequence:
            return ''.join(sequence)
        else:
            logging.warning(f"No sequence found in {filename} for chain {chain}.")
            return None
    except FileNotFoundError:
        logging.error(f"File not found: {filename}")
        return None
    except Exception as e:
        logging.error(f"Error processing {filename}: {e}")
        return None

def process_pdb_entry(entry: List[str], sc_directory: Path, chain: str = "A") -> Optional[List[str]]:
    """
    Process a single pdb entry to extract sequence and size.

    Args:
        entry (List[str]): A list containing the required fields from the .sc file.
        sc_directory (Path): Directory containing PDB files.
        chain (str): Chain identifier to extract the sequence from.

    Returns:
        Optional[List[str]]: A list containing all original fields plus Sequence and Size,
                             or None if processing fails.
    """
    if len(entry) < 5:
        logging.warning("Invalid entry with insufficient columns.")
        return None
    binder_aligned_rmsd, pae_binder, pae_interaction, plddt_binder, description = entry
    pdb_id = description  # Assuming 'description' is the PDB ID without '.pdb'
    pdb_filename = sc_directory / f"{pdb_id}.pdb"
    sequence = extract_sequence(pdb_filename, chain=chain)
    if sequence:
        size = len(sequence)
        return [binder_aligned_rmsd, pae_binder, pae_interaction, plddt_binder, description, sequence, str(size)]
    else:
        logging.info(f"No sequence extracted for PDB ID: {pdb_id}")
        return None

def process_pdb_entries_parallel(filtered_entries: List[List[str]],
                                 output_file: Path,
                                 sc_directory: Path,
                                 chain: str,
                                 max_workers: int) -> None:
    """
    Process filtered pdb entries in parallel to extract sequences and generate binder.csv.

    Args:
        filtered_entries (List[List[str]]): Filtered data from .sc files.
        output_file (Path): Path to the output binder.csv file.
        sc_directory (Path): Directory containing PDB files.
        chain (str): Chain identifier to extract sequences from.
        max_workers (int): The maximum number of threads to use for parallel processing.
    """
    # Define CSV headers
    csv_headers = ["binder_aligned_rmsd", "pae_binder", "pae_interaction", "plddt_binder", "description", "sequence", "size"]

    # Write CSV header
    with output_file.open("w", newline='') as csvfile:
        writer = csv.writer(csvfile)
        writer.writerow(csv_headers)

    logging.info(f"Starting processing with {max_workers} workers...")

    with ThreadPoolExecutor(max_workers=max_workers) as executor:
        # Submit all tasks
        futures = {
            executor.submit(process_pdb_entry, entry, sc_directory, chain): entry
            for entry in filtered_entries
        }

        processed_count = 0
        total_entries = len(futures)
        for future in as_completed(futures):
            result = future.result()
            if result:
                with output_file.open("a", newline='') as csvfile:
                    writer = csv.writer(csvfile)
                    writer.writerow(result)
                processed_count += 1
                if processed_count % 100 == 0 or processed_count == total_entries:
                    logging.info(f"Processed {processed_count}/{total_entries} entries.")

    logging.info(f"Processing completed successfully. {processed_count}/{total_entries} entries written to {output_file}.")

def main():
    # Set up argument parser
    parser = argparse.ArgumentParser(description="Process PDB files to generate binder.csv without intermediate files.")
    parser.add_argument('--settings_file', required=True, help="Path to compute_settings.txt file")
    parser.add_argument('--output_csv', type=str, default="binders.csv", help='Path to the output binders.csv file.')
    parser.add_argument('--chain', type=str, default="A", help='Chain identifier to extract sequences from.')
    parser.add_argument('--sc_directory', type=str, default=".", help='Directory containing out_*.sc files and PDB files.')
    
    args = parser.parse_args()

    settings_file = Path(args.settings_file)
    try:
        settings = parse_settings_file(settings_file)
    except Exception as e:
        logging.error(f"Failed to parse settings file: {e}")
        sys.exit(1)

    # Extract max_workers from settings
    try:
        max_workers = int(settings.get('max_workers', 48))  # Default to 48 if not specified
        logging.info(f"Using max_workers: {max_workers}")
    except ValueError:
        logging.error("Invalid value for max_workers in settings file. It must be an integer.")
        sys.exit(1)

    # Load Python module
    load_python_module()

    sc_directory = Path(args.sc_directory)
    if not sc_directory.is_dir():
        logging.error(f"Specified sc_directory does not exist or is not a directory: {sc_directory}")
        sys.exit(1)

    # Combine out_*.sc files into memory
    combined_lines = combine_sc_files(sc_directory)

    # Filter lines based on pae_interaction
    filtered_entries = filter_sc_lines(combined_lines)

    if not filtered_entries:
        logging.error("No entries found after filtering. Exiting.")
        sys.exit(1)

    # Process filtered entries in parallel and generate binder.csv
    output_csv = Path(args.output_csv)
    process_pdb_entries_parallel(
        filtered_entries=filtered_entries,
        output_file=output_csv,
        sc_directory=sc_directory,
        chain=args.chain,
        max_workers=max_workers
    )

    logging.info("All tasks completed successfully.")

if __name__ == "__main__":
    main()


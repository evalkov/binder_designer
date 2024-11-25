#!/usr/bin/env python3

import os
import sys
import csv
import argparse
import logging
from pathlib import Path
import xml.etree.ElementTree as ET
from concurrent.futures import ProcessPoolExecutor, as_completed
from typing import List, Optional, Dict
import pandas as pd

# Constants
DEFAULT_LOG_FILE = "process_salt_bridges_and_filter_sc.log"
DEFAULT_OUTPUT_CSV = "qualifying_salt_bridges.csv"
DEFAULT_FILTERED_SC_FILE = "filtered_out.sc"
DEFAULT_SETTINGS_FILE = "compute_settings.txt"


def setup_logging(log_file: str = DEFAULT_LOG_FILE) -> None:
    """
    Configure logging to output messages to both a log file and the console.

    Args:
        log_file (str): Path to the log file.
    """
    logging.basicConfig(
        level=logging.INFO,
        format='%(asctime)s [%(levelname)s] %(message)s',
        handlers=[
            logging.FileHandler(log_file, mode='w'),
            logging.StreamHandler(sys.stdout)
        ]
    )


def parse_arguments() -> argparse.Namespace:
    """
    Parse command-line arguments.

    Returns:
        argparse.Namespace: Parsed command-line arguments.
    """
    parser = argparse.ArgumentParser(
        description="Process PISA XML files to extract qualifying salt bridges and filter a .sc file based on the results."
    )
    
    # Arguments for XML processing
    parser.add_argument(
        '--xml_directory',
        required=True,
        type=str,
        help="Path to the directory containing PISA XML files."
    )
    parser.add_argument(
        '--chain_a',
        type=str,
        default='A',
        help="Chain ID for Chain A (default: A)."
    )
    parser.add_argument(
        '--chain_b',
        type=str,
        default='B',
        help="Chain ID for Chain B (default: B)."
    )
    parser.add_argument(
        '--residue_offset',
        required=True,
        type=int,
        help="Residue offset to add to the first residue of Chain B."
    )
    
    # Arguments for .sc filtering
    parser.add_argument(
        '--af2_score',
        required=True,
        metavar='AF2_SCORE_FILE',
        help='Path to the input .sc file (e.g., out.sc).'
    )
    parser.add_argument(
        '--filtered_sc_output',
        type=str,
        default=DEFAULT_FILTERED_SC_FILE,
        metavar='FILTERED_SC_OUTPUT',
        help='Path to the output filtered .sc file. Default is "filtered_out.sc".'
    )
    
    # Common arguments
    parser.add_argument(
        '--output_csv',
        type=str,
        default=DEFAULT_OUTPUT_CSV,
        metavar='OUTPUT_CSV',
        help=f"Path to the output CSV file (default: {DEFAULT_OUTPUT_CSV})."
    )
    parser.add_argument(
        '--settings_file',
        type=str,
        default=DEFAULT_SETTINGS_FILE,
        help=f"Path to the compute settings file (default: {DEFAULT_SETTINGS_FILE})."
    )
    parser.add_argument(
        '--max_workers',
        type=int,
        default=4,
        help="Number of parallel workers (processes) to use (default: 4)."
    )
    
    return parser.parse_args()


def read_compute_settings(settings_file: Path) -> Dict[str, Optional[int]]:
    """
    Reads compute settings from a file.

    Parameters:
        settings_file (Path): Path to the settings file.

    Returns:
        dict: A dictionary containing compute settings (e.g., max_workers, batch_size).
    """
    settings = {'max_workers': 4, 'batch_size': None}  # Defaults
    if settings_file.is_file():
        with settings_file.open('r') as f:
            lines = f.readlines()
        # Only parse lines under "Recommended settings:"
        start_parsing = False
        for line in lines:
            if "Recommended settings:" in line:
                start_parsing = True
                continue
            if start_parsing:
                if 'max_workers' in line:
                    try:
                        settings['max_workers'] = int(line.split(':')[1].strip())
                        logging.info(f"Setting max_workers to {settings['max_workers']} from settings file.")
                    except ValueError:
                        logging.warning("Invalid value for max_workers in settings file. Using default.")
                if 'batch_size' in line:
                    try:
                        settings['batch_size'] = int(line.split(':')[1].strip())
                        logging.info(f"Setting batch_size to {settings['batch_size']} from settings file.")
                    except ValueError:
                        logging.warning("Invalid value for batch_size in settings file. Using default.")
    else:
        logging.warning(f"Settings file {settings_file} does not exist. Using default settings.")
    return settings


def process_single_file(filepath: Path, target_chain_A: str, target_chain_B: str, residue_offset: int) -> List[List[str]]:
    """
    Processes a single XML file and extracts salt bridges.

    Parameters:
        filepath (Path): Path to the XML file.
        target_chain_A (str): Chain ID for Chain A.
        target_chain_B (str): Chain ID for Chain B.
        residue_offset (int): Offset to add to the first residue of Chain B.

    Returns:
        list: List of salt bridge data in the format:
              [formatted_pdb_id, Residue1, Residue2, Distance]
    """
    results = []
    filename = filepath.name
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
            logging.warning(f"Chain {target_chain_B} not found in {filename}. Skipping interface.")
            continue

        first_seq_num_B = min(chain_residues[target_chain_B])
        target_residue_B = first_seq_num_B + residue_offset

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

            # Check if the bond matches the target criteria
            condition1 = (chain1_elem.text.strip() == target_chain_B and
                          seqnum1 == target_residue_B and
                          chain2_elem.text.strip() == target_chain_A)

            condition2 = (chain2_elem.text.strip() == target_chain_B and
                          seqnum2 == target_residue_B and
                          chain1_elem.text.strip() == target_chain_A)

            if condition1 or condition2:
                res1 = f"{chain1_elem.text.strip()}/{res1_elem.text.strip()}{seqnum1}"
                res2 = f"{chain2_elem.text.strip()}/{res2_elem.text.strip()}{seqnum2}"
                dist_rounded = round(dist, 2)

                formatted_pdb_id = pdb_id  # Already lowercase
                logging.info(f"{formatted_pdb_id}: {res1} ↔ {res2} at {dist_rounded} Å")

                results.append([formatted_pdb_id, res1, res2, dist_rounded])

    return results


def find_pdb_files_with_specific_salt_bridge(xml_directory: Path,
                                             target_chain_A: str,
                                             target_chain_B: str,
                                             output_csv: Path,
                                             max_workers: int,
                                             residue_offset: int) -> None:
    """
    Processes XML files in parallel to extract salt bridges.

    Parameters:
        xml_directory (Path): Path to the directory containing XML files.
        target_chain_A (str): Chain ID for Chain A.
        target_chain_B (str): Chain ID for Chain B.
        output_csv (Path): Path to the output CSV file.
        max_workers (int): Number of parallel workers (processes).
        residue_offset (int): Offset to add to the first residue of Chain B.
    """
    results = []

    # List all XML files in the directory
    filepaths = [f for f in xml_directory.glob('*.xml') if f.is_file()]
    logging.info(f"Found {len(filepaths)} XML files in {xml_directory}.")

    if not filepaths:
        logging.warning(f"No XML files found in {xml_directory}. Exiting XML processing step.")
        return

    with ProcessPoolExecutor(max_workers=max_workers) as executor:
        # Submit all tasks
        future_to_file = {
            executor.submit(process_single_file, filepath, target_chain_A, target_chain_B, residue_offset): filepath
            for filepath in filepaths
        }

        for future in as_completed(future_to_file):
            filepath = future_to_file[future]
            try:
                file_results = future.result()
                if file_results:
                    results.extend(file_results)
            except Exception as e:
                logging.error(f"Error processing file {filepath.name}: {e}")

    # Write results to CSV
    if results:
        try:
            with output_csv.open("w", newline='') as csvfile:
                csv_writer = csv.writer(csvfile)
                csv_writer.writerow(['PDB_ID', 'Residue1', 'Residue2', 'Distance (Å)'])
                csv_writer.writerows(results)
            logging.info(f"Salt bridges successfully written to {output_csv}")
        except Exception as e:
            logging.error(f"Failed to write to {output_csv}: {e}")
    else:
        logging.info("No qualifying salt bridges found. CSV file not created.")


def filter_sc_file(salt_bridges_csv: Path, sc_file: Path, filtered_sc_output: Path) -> None:
    """
    Filters the .sc file based on PDB_IDs from the salt bridges CSV.

    Parameters:
        salt_bridges_csv (Path): Path to the qualifying_salt_bridges.csv file.
        sc_file (Path): Path to the input .sc file.
        filtered_sc_output (Path): Path to the output filtered .sc file.
    """
    # Step 1: Load PDB_IDs from qualifying_salt_bridges.csv
    logging.info(f"Reading {salt_bridges_csv}...")
    if not salt_bridges_csv.is_file():
        logging.error(f"Salt bridges CSV file does not exist: {salt_bridges_csv}")
        sys.exit(1)

    try:
        salt_bridges_df = pd.read_csv(salt_bridges_csv)
    except Exception as e:
        logging.error(f"Failed to read {salt_bridges_csv}: {e}")
        sys.exit(1)

    # Check if 'PDB_ID' column exists
    if 'PDB_ID' not in salt_bridges_df.columns:
        logging.error(f"'PDB_ID' column not found in {salt_bridges_csv}.")
        sys.exit(1)

    # Extract unique PDB_IDs and ensure they are clean
    pdb_ids = set(salt_bridges_df["PDB_ID"].astype(str).str.strip())
    logging.info(f"Extracted {len(pdb_ids)} unique PDB_IDs from CSV.")

    # Step 2: Read and process the .sc file
    logging.info(f"Reading {sc_file}...")
    if not sc_file.is_file():
        logging.error(f".sc file does not exist: {sc_file}")
        sys.exit(1)

    try:
        with sc_file.open("r") as f:
            lines = f.readlines()
    except Exception as e:
        logging.error(f"Failed to read {sc_file}: {e}")
        sys.exit(1)

    # Separate header and data lines
    if not lines:
        logging.error(f".sc file {sc_file} is empty.")
        sys.exit(1)

    header_line = lines[0]  # The first line is the header
    data_lines = lines[1:]  # The remaining lines are data

    logging.info(f"Found 1 header line and {len(data_lines)} data lines in {sc_file}.")

    # Step 3: Filter data_lines based on matching description field
    logging.info("Filtering .sc data lines based on PDB_IDs...")
    filtered_lines = []
    for line in data_lines:
        # Split line into fields and get the last field as description
        parts = line.strip().split()
        if len(parts) > 0:
            description = parts[-1].strip()  # Last field is the description
            if description in pdb_ids:
                filtered_lines.append(line)  # Add line if description matches PDB_ID

    logging.info(f"Filtered {len(filtered_lines)} matching lines from .sc file.")

    # Step 4: Write the filtered output
    logging.info(f"Writing filtered .sc output to {filtered_sc_output}...")
    try:
        with filtered_sc_output.open("w") as f:
            f.write(header_line)  # Write the header
            f.writelines(filtered_lines)  # Write filtered data
        logging.info(f"Filtered .sc file successfully saved to {filtered_sc_output}")
    except Exception as e:
        logging.error(f"Failed to write to {filtered_sc_output}: {e}")
        sys.exit(1)


def main():
    # Parse command-line arguments
    args = parse_arguments()

    # Initialize logging
    setup_logging()

    # Convert arguments to Path objects
    xml_directory = Path(args.xml_directory)
    output_csv = Path(args.output_csv)
    filtered_sc_output = Path(args.filtered_sc_output)
    settings_file = Path(args.settings_file)
    af2_score_file = Path(args.af2_score)
    residue_offset = args.residue_offset

    # Validate xml_directory
    if not xml_directory.is_dir():
        logging.critical(f"The directory {xml_directory} does not exist.")
        print(f"Error: The directory {xml_directory} does not exist. Check the path and try again.")
        sys.exit(1)

    # Read compute settings
    compute_settings = read_compute_settings(settings_file)
    max_workers = compute_settings.get('max_workers', args.max_workers)

    # Log the number of workers being used
    logging.info(f"Using max_workers: {max_workers}")

    # Step 1: Process XML files to extract salt bridges and generate CSV
    logging.info("Starting salt bridge extraction from XML files...")
    find_pdb_files_with_specific_salt_bridge(
        xml_directory=xml_directory,
        target_chain_A=args.chain_a,
        target_chain_B=args.chain_b,
        output_csv=output_csv,
        max_workers=max_workers,
        residue_offset=residue_offset
    )
    logging.info("Salt bridge extraction completed.")

    # Check if the CSV was created and has entries
    if not output_csv.is_file():
        logging.error(f"CSV file {output_csv} was not created. Skipping .sc filtering step.")
        sys.exit(1)

    # Step 2: Filter .sc file based on the generated CSV
    logging.info("Starting .sc file filtering based on extracted salt bridges...")
    filter_sc_file(
        salt_bridges_csv=output_csv,
        sc_file=af2_score_file,
        filtered_sc_output=filtered_sc_output
    )
    logging.info("All processes completed successfully.")


if __name__ == "__main__":
    main()


#!/usr/bin/env python3

import os
import shutil
import subprocess
import xml.etree.ElementTree as ET
import argparse
import logging
from pathlib import Path
import sys
from concurrent.futures import ThreadPoolExecutor, as_completed
import threading

def setup_logging():
    """
    Set up logging to output to both console and a log file.
    """
    log_file = 'pisa_analysis.log'  # Default log file name

    logging.basicConfig(
        level=logging.INFO,  # Set to DEBUG for more detailed logs
        format='%(asctime)s [%(levelname)s] %(message)s',
        handlers=[
            logging.FileHandler(log_file),
            logging.StreamHandler(sys.stdout)
        ]
    )
    logging.info("Logging is set up.")

def parse_arguments():
    """
    Parse command-line arguments.

    Returns:
        argparse.Namespace: Parsed arguments.
    """
    parser = argparse.ArgumentParser(description="Process PDB files with PISA analysis.")
    parser.add_argument('--settings_file', '-c', required=True, help="Path to compute_settings.txt file")
    return parser.parse_args()

def load_ccp4_environment():
    """
    Source the CCP4 environment persistently.
    """
    try:
        logging.info("Loading CCP4 environment...")
        result = subprocess.run(
            "module load ccp4/8.0.010 && env",
            shell=True,
            check=True,
            stdout=subprocess.PIPE,
            executable="/bin/bash"
        )
        env_vars = {}
        for line in result.stdout.decode("utf-8").splitlines():
            if '=' in line:
                key, value = line.split('=', 1)
                env_vars[key] = value

        os.environ.update(env_vars)
        logging.info("CCP4 environment variables loaded successfully.")
    except subprocess.CalledProcessError as e:
        logging.error("Failed to load CCP4 environment. Ensure the module system is available and the CCP4 module is valid.")
        raise RuntimeError("Failed to load CCP4 environment.") from e

    os.environ['CCP4'] = '/mnt/nasapps/production/ccp4/8.0.010'
    ccp4_bin = os.path.join(os.environ['CCP4'], 'bin')
    os.environ['PATH'] = f"{ccp4_bin}:{os.environ.get('PATH', '')}"
    logging.info(f"Updated PATH with CCP4 binaries: {ccp4_bin}")

    if shutil.which("pisa") is None:
        logging.error("PISA executable not found. Ensure CCP4 is loaded correctly and added to the PATH.")
        raise RuntimeError("PISA executable not found.")
    else:
        logging.info("PISA executable found.")

def read_settings(settings_file):
    """
    Read and parse the settings file to extract configuration parameters.

    Args:
        settings_file (str): Path to the settings file.

    Returns:
        dict: Dictionary containing configuration parameters.
    """
    if not os.path.isfile(settings_file):
        logging.error(f"The specified settings file {settings_file} does not exist.")
        raise ValueError(f"The specified settings file {settings_file} does not exist.")
    
    settings = {}
    try:
        with open(settings_file, 'r') as f:
            lines = f.readlines()
        for line in lines:
            line = line.strip()
            if not line or line.startswith('#'):
                continue  # Skip empty lines and comments
            if ':' in line:
                key, value = line.split(':', 1)
                settings[key.strip()] = value.strip()
        logging.info(f"Settings loaded from {settings_file}: {settings}")
    except Exception as e:
        logging.error(f"Failed to read settings file {settings_file}: {e}")
        raise

    # Validate required settings
    if 'max_workers' not in settings:
        logging.error(f"'max_workers' not found in settings file {settings_file}.")
        raise ValueError(f"'max_workers' not found in settings file {settings_file}.")

    try:
        settings['max_workers'] = int(settings['max_workers'])
    except ValueError:
        logging.error(f"'max_workers' must be an integer in settings file {settings_file}.")
        raise ValueError(f"'max_workers' must be an integer in settings file {settings_file}.")

    if 'batch_size' in settings:
        try:
            settings['batch_size'] = int(settings['batch_size'])
        except ValueError:
            logging.error(f"'batch_size' must be an integer in settings file {settings_file}.")
            raise ValueError(f"'batch_size' must be an integer in settings file {settings_file}.")
    else:
        settings['batch_size'] = None  # or set a default value if desired

    return settings

def initialize_parameters():
    """
    Initialize residue_counter and chain_id with default values if they are None.

    Returns:
        tuple: (residue_counter, chain_id)
    """
    residue_counter = None
    chain_id = None

    if residue_counter is None:
        residue_counter = 0
    if chain_id is None:
        chain_id = 'A'

    if residue_counter == 0:
        logging.info(f"Analyzing chain {chain_id} with no residue restriction")
    else:
        logging.info(f"Analyzing chain {chain_id} and only after residue {residue_counter}")

    return residue_counter, chain_id

def parse_xml(file, chain_id, residue_counter):
    """
    Parse the XML file generated by PISA to extract residues and bond counts.

    Args:
        file (str): Path to the XML file.
        chain_id (str): Chain identifier to analyze.
        residue_counter (int): Residue counter for analysis.

    Returns:
        tuple: (residues list, h_bonds_count, salt_bridges_count)
    """
    try:
        tree = ET.parse(file)
        root = tree.getroot()
        logging.info(f"Successfully parsed XML file: {file}")
    except ET.ParseError as e:
        logging.error(f"Error parsing {file}: {e}")
        return [], 0, 0
    except Exception as e:
        logging.error(f"Unexpected error parsing {file}: {e}")
        return [], 0, 0

    residues = []
    h_bonds_count = 0
    salt_bridges_count = 0

    for molecule in root.findall('.//molecule'):
        current_chain_id = molecule.find('chain_id').text
        for residue in molecule.findall('.//residue'):
            res_name = residue.find('name').text
            seq_num = int(residue.find('seq_num').text)
            asa = float(residue.find('asa').text)
            bsa = float(residue.find('bsa').text)
            solv_en = float(residue.find('solv_en').text)
            residues.append((current_chain_id, res_name, seq_num, asa, bsa, solv_en))

    for interface in root.findall('.//interface'):
        h_bonds = interface.find('h-bonds')
        for bond in h_bonds.findall('bond'):
            if bond.find('chain-1').text == chain_id or bond.find('chain-2').text == chain_id:
                seqnum_1 = int(bond.find('seqnum-1').text)
                seqnum_2 = int(bond.find('seqnum-2').text)
                if bond.find('chain-1').text == chain_id and seqnum_1 >= residue_counter:
                    h_bonds_count += 1
                if bond.find('chain-2').text == chain_id and seqnum_2 >= residue_counter:
                    h_bonds_count += 1

        salt_bridges = interface.find('salt-bridges')
        for bond in salt_bridges.findall('bond'):
            if bond.find('chain-1').text == chain_id or bond.find('chain-2').text == chain_id:
                seqnum_1 = int(bond.find('seqnum-1').text)
                seqnum_2 = int(bond.find('seqnum-2').text)
                if bond.find('chain-1').text == chain_id and seqnum_1 >= residue_counter:
                    salt_bridges_count += 1
                if bond.find('chain-2').text == chain_id and seqnum_2 >= residue_counter:
                    salt_bridges_count += 1

    return residues, h_bonds_count, salt_bridges_count

def calculate_buried_area_percentage(residues):
    """
    Calculate the buried area percentage for each residue.

    Args:
        residues (list): List of residues.

    Returns:
        list: List of tuples with buried area data.
    """
    output = []
    for residue in residues:
        current_chain_id, res_name, seq_num, asa, bsa, solv_en = residue
        if bsa != 0 and asa != 0:
            percentage = (bsa / asa) * 100
            bars = int(percentage / 10)
            bar_str = "|" * bars
            output.append((current_chain_id, res_name, seq_num, asa, bsa, solv_en, percentage, bar_str))
    return output

def calculate_total_bsa(residues, chain_id, residue_counter):
    """
    Calculate the total buried surface area (BSA) score.

    Args:
        residues (list): List of residues.
        chain_id (str): Chain identifier to analyze.
        residue_counter (int): Residue counter for analysis.

    Returns:
        int: Total BSA score.
    """
    total_bars = 0
    for residue in residues:
        current_chain_id, _, seq_num, asa, bsa, _ = residue
        if current_chain_id == chain_id and seq_num >= residue_counter and asa != 0:
            percentage = (bsa / asa) * 100
            bars = int(percentage / 10)
            total_bars += bars
    return total_bars

def process_pdb_file(file, chain_id, residue_counter):
    """
    Process a single PDB file with PISA analysis.

    Args:
        file (str): Path to the PDB file.
        chain_id (str): Chain identifier to analyze.
        residue_counter (int): Residue counter for analysis.

    Returns:
        tuple: (base_filename, total_bsa_score, h_bonds_count, salt_bridges_count)
    """
    base_filename = os.path.splitext(file)[0]
    try:
        logging.info(f"Starting PISA analysis for {file}")
        
        # Run PISA analyze
        subprocess.run(['pisa', file, '-analyse', file], stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL, check=True)
        logging.info(f"PISA analysis completed for {file}")
        
        # Run PISA XML generation
        xml_output = f"{base_filename}.xml"
        with open(xml_output, 'w') as xml_file:
            subprocess.run(['pisa', file, '-xml', 'interfaces'], stdout=xml_file, stderr=subprocess.DEVNULL, check=True)
        logging.info(f"PISA XML generated for {file} as {xml_output}")
        
        # Run PISA erase
        subprocess.run(['pisa', file, '-erase'], stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL, check=True)
        logging.info(f"PISA erase completed for {file}")
        
    except subprocess.CalledProcessError as e:
        logging.error(f"PISA analysis failed for {file}: {e}")
        return base_filename, 0, 0, 0

    xml_file = f"{base_filename}.xml"
    residues, h_bonds_count, salt_bridges_count = parse_xml(xml_file, chain_id, residue_counter)
    buried_area_percentage_data = calculate_buried_area_percentage(residues)
    total_bsa_score = calculate_total_bsa(residues, chain_id, residue_counter)

    return base_filename, total_bsa_score, h_bonds_count, salt_bridges_count

def main():
    # Set up logging
    setup_logging()
    logging.info("=== Starting PISA analysis process ===")

    # Parse command-line arguments
    args = parse_arguments()
    settings_file = args.settings_file

    # Load CCP4 environment
    load_ccp4_environment()

    # Read settings
    settings = read_settings(settings_file)
    max_workers = settings['max_workers']
    batch_size = settings.get('batch_size', None)  # Currently unused

    # Initialize parameters
    residue_counter, chain_id = initialize_parameters()

    # Get list of PDB files
    pdb_dir = os.getcwd()
    pdb_files = [file for file in os.listdir(pdb_dir) if file.endswith('.pdb')]
    total_files = len(pdb_files)
    logging.info(f"Found {total_files} .pdb files in '{pdb_dir}'.")

    if total_files == 0:
        logging.warning(f"No .pdb files found in the directory '{pdb_dir}'. Exiting.")
        sys.exit(0)

    logging.info(f"Analyzing chain {chain_id} with residue counter set to {residue_counter}.")

    results = []
    processed_count = 0
    last_reported_percentage = 0

    # Lock for thread-safe updates
    lock = threading.Lock()

    with ThreadPoolExecutor(max_workers=max_workers) as executor:
        futures = {executor.submit(process_pdb_file, file, chain_id, residue_counter): file for file in pdb_files}
        logging.info(f"Submitted {len(futures)} PDB files for processing with max_workers={max_workers}.")

        for idx, future in enumerate(as_completed(futures), start=1):
            file = futures[future]
            try:
                result = future.result()
                results.append(result)
                processed_count += 1
                percentage = (processed_count / total_files) * 100
                if int(percentage) // 10 > last_reported_percentage // 10:
                    last_reported_percentage = int(percentage)
                    logging.info(f"Processed {processed_count}/{total_files} files ({last_reported_percentage}% complete)")
            except Exception as e:
                logging.error(f"Error processing {file}: {e}")

    # Move XML files to 'pisa_xml_files' directory
    xml_dir = 'pisa_xml_files'
    try:
        os.makedirs(xml_dir, exist_ok=True)
        logging.info(f"Created/made sure the directory '{xml_dir}' exists for XML files.")
    except Exception as e:
        logging.error(f"Failed to create directory '{xml_dir}': {e}")
        sys.exit(1)

    moved_files = 0
    for file in os.listdir('.'):
        if file.endswith('.xml'):
            try:
                shutil.move(file, os.path.join(xml_dir, file))
                moved_files += 1
                logging.info(f"Moved XML file '{file}' to '{xml_dir}'.")
            except Exception as e:
                logging.error(f"Failed to move XML file '{file}' to '{xml_dir}': {e}")

    logging.info(f"Moved {moved_files} XML files to '{xml_dir}'.")

    # Write results to 'contacts.csv'
    contacts_csv = 'contacts.csv'
    try:
        with open(contacts_csv, 'w') as f:
            f.write("binder,bsa_score,salt_bridges,h_bonds\n")
            for result in results:
                base_filename, total_bsa_score, salt_bridges_count, h_bonds_count = result
                f.write(f"{base_filename},{total_bsa_score},{salt_bridges_count},{h_bonds_count}\n")
        logging.info(f"Successfully wrote results to '{contacts_csv}'.")
    except Exception as e:
        logging.error(f"Failed to write results to '{contacts_csv}': {e}")
        sys.exit(1)

    logging.info(f"Processing complete. {len(results)} files out of {total_files} were successfully processed.")
    logging.info("=== PISA analysis process completed successfully ===")

if __name__ == "__main__":
    main()


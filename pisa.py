import os
import shutil
import subprocess
import xml.etree.ElementTree as ET
from concurrent.futures import ThreadPoolExecutor, as_completed

# Read compute settings from the file
settings_file = 'compute_settings.txt'
with open(settings_file, 'r') as f:
    lines = f.readlines()

# Extract max_workers and batch_size from the settings file
for line in lines:
    if 'max_workers' in line:
        max_workers = int(line.split(':')[1].strip())
    if 'batch_size' in line:
        batch_size = int(line.split(':')[1].strip())

# Set the residue counter and chain ID with default values
residue_counter = None
chain_id = None

# Check if residue_counter and chain_id are specified
if residue_counter is None:
    residue_counter = 0  # No residue restriction
if chain_id is None:
    chain_id = 'A'  # Default chain ID

# Print the chain and residue information
if residue_counter == 0:
    print(f"Analyzing chain {chain_id} with no residue restriction")
else:
    print(f"Analyzing chain {chain_id} and only after residue {residue_counter}")

# Directory containing PDB files
pdb_dir = '.'
pdb_files = [file for file in os.listdir(pdb_dir) if file.endswith('.pdb')]

# Function to parse XML and extract relevant data
def parse_xml(file, chain_id, residue_counter):
    try:
        tree = ET.parse(file)
        root = tree.getroot()
    except ET.ParseError as e:
        print(f"Error parsing {file}: {e}")
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
        # Hydrogen bonds
        h_bonds = interface.find('h-bonds')
        for bond in h_bonds.findall('bond'):
            if bond.find('chain-1').text == chain_id or bond.find('chain-2').text == chain_id:
                seqnum_1 = int(bond.find('seqnum-1').text)
                seqnum_2 = int(bond.find('seqnum-2').text)
                if bond.find('chain-1').text == chain_id and seqnum_1 >= residue_counter:
                    h_bonds_count += 1
                elif bond.find('chain-2').text == chain_id and seqnum_2 >= residue_counter:
                    h_bonds_count += 1

        # Salt bridges
        salt_bridges = interface.find('salt-bridges')
        for bond in salt_bridges.findall('bond'):
            if bond.find('chain-1').text == chain_id or bond.find('chain-2').text == chain_id:
                seqnum_1 = int(bond.find('seqnum-1').text)
                seqnum_2 = int(bond.find('seqnum-2').text)
                if bond.find('chain-1').text == chain_id and seqnum_1 >= residue_counter:
                    salt_bridges_count += 1
                elif bond.find('chain-2').text == chain_id and seqnum_2 >= residue_counter:
                    salt_bridges_count += 1

    return residues, h_bonds_count, salt_bridges_count

# Function to calculate buried area percentage
def calculate_buried_area_percentage(residues):
    output = []
    for residue in residues:
        current_chain_id, res_name, seq_num, asa, bsa, solv_en = residue
        if bsa != 0 and asa != 0:
            percentage = (bsa / asa) * 100
            bars = int(percentage / 10)
            bar_str = "|" * bars
            output.append((current_chain_id, res_name, seq_num, asa, bsa, solv_en, percentage, bar_str))
    return output

# Function to calculate the total BSA score for residues in the specified chain with sequence number >= residue_counter
def calculate_total_bsa(residues, chain_id, residue_counter):
    total_bars = 0
    for residue in residues:
        current_chain_id, _, seq_num, asa, bsa, _ = residue
        if current_chain_id == chain_id and seq_num >= residue_counter and asa != 0:
            percentage = (bsa / asa) * 100
            bars = int(percentage / 10)
            total_bars += bars
    return total_bars

# Function to process a single PDB file
def process_pdb_file(file, chain_id, residue_counter):
    base_filename = os.path.splitext(file)[0]

    # Run PISA analysis and create an XML output file
    subprocess.run(['pisa', file, '-analyse', file], stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)
    subprocess.run(['pisa', file, '-xml', 'interfaces'], stdout=open(f'{base_filename}.xml', 'w'), stderr=subprocess.DEVNULL)
    subprocess.run(['pisa', file, '-erase'], stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)

    # Load and process the XML file
    xml_file = f'{base_filename}.xml'
    residues, h_bonds_count, salt_bridges_count = parse_xml(xml_file, chain_id, residue_counter)

    # Calculate buried area percentage
    buried_area_percentage_data = calculate_buried_area_percentage(residues)

    # Calculate total BSA score for the specified chain
    total_bsa_score = calculate_total_bsa(residues, chain_id, residue_counter)

    return base_filename, total_bsa_score, h_bonds_count, salt_bridges_count

# Process files in batches
results = []
with ThreadPoolExecutor(max_workers=max_workers) as executor:
    futures = {executor.submit(process_pdb_file, file, chain_id, residue_counter): file for file in pdb_files}

    for future in as_completed(futures):
        file = futures[future]
        try:
            result = future.result()
            results.append(result)
        except Exception as e:
            print(f"Error processing {file}: {e}")

# Write all results to the output file
with open('contacts.csv', 'w') as f:
    f.write("binder,bsa_score,salt_bridges,h_bonds\n")
    for result in results:
        base_filename, total_bsa_score, salt_bridges_count, h_bonds_count = result
        f.write(f"{base_filename},{total_bsa_score},{salt_bridges_count},{h_bonds_count}\n")

# Create a directory for XML files if it doesn't exist
xml_dir = 'pisa_xml_files'
os.makedirs(xml_dir, exist_ok=True)

# Move all XML files to the new directory
for file in os.listdir('.'):
    if file.endswith('.xml'):
        shutil.move(file, os.path.join(xml_dir, file))

print(f"Processed {len(results)} files.")

import os
import csv
import xml.etree.ElementTree as ET
from concurrent.futures import ThreadPoolExecutor, as_completed
import psutil

CCP4 = "/mnt/nasapps/production/ccp4/8.0.010"
INPUT_DIR = "/mnt/beegfs/valkov/rfdiffusion/binder781-model1/binder781-model1-TFldTo5K"
OUTPUT_DIR = "/scratch/cluster_scratch/valkove2/top_scoring_solutions"
residue_counter = 44
batch_size = 2054  # Adjust based on memory and file size
max_workers = 48  # Adjust based on available cpu cores

os.chdir(INPUT_DIR)

header = ["binder", "bsa_score", "salt_bridges", "h_bonds"]
contacts_file = os.path.join(OUTPUT_DIR, "contacts.csv")

def is_well_formed(xml_file):
    try:
        ET.parse(xml_file)
        return True
    except ET.ParseError:
        return False

def process_file(file):
    base_filename = os.path.splitext(file)[0]
    xml_file = os.path.join(OUTPUT_DIR, f"{base_filename}.xml")
    
    # Run PISA commands
    os.system(f"pisa {file} -analyse {file} > /dev/null 2>&1")
    os.system(f"pisa {file} -xml interfaces > {xml_file}")
    os.system(f"pisa {file} -erase > /dev/null 2>&1")
    
    # Check if the XML file is well-formed before parsing
    if not is_well_formed(xml_file):
        print(f"Error: {xml_file} is not well-formed. Skipping this file.")
        return None
    
    total_bsa = 0
    total_salt_bridges = 0
    total_hbonds = 0

    with open(os.path.join(OUTPUT_DIR, f"{base_filename}_sa.csv"), mode='w', newline='') as sa_file, \
         open(os.path.join(OUTPUT_DIR, f"{base_filename}_salt_bridges.csv"), mode='w', newline='') as sb_file, \
         open(os.path.join(OUTPUT_DIR, f"{base_filename}_hbonds.csv"), mode='w', newline='') as hb_file:
        
        sa_writer = csv.writer(sa_file)
        sb_writer = csv.writer(sb_file)
        hb_writer = csv.writer(hb_file)
        
        sa_writer.writerow(["Chain", "Residue Name", "Sequence Number", "ASA", "BSA", "Solvent Energy"])
        sb_writer.writerow(["Chain 1", "Residue 1", "Sequence Number 1", "Atom Name 1", "Chain 2", "Residue 2", "Sequence Number 2", "Atom Name 2", "Distance"])
        hb_writer.writerow(["Chain 1", "Residue 1", "Sequence Number 1", "Atom Name 1", "Chain 2", "Residue 2", "Sequence Number 2", "Atom Name 2", "Distance"])
        
        context = ET.iterparse(xml_file, events=("start", "end"))
        context = iter(context)
        event, root = next(context)
        
        for event, elem in context:
            if event == "end" and elem.tag == "molecule":
                chain_id = elem.find('chain_id').text
                for residue in elem.findall('.//residue'):
                    res_name = residue.find('name').text
                    seq_num = int(residue.find('seq_num').text)
                    asa = float(residue.find('asa').text)
                    bsa = float(residue.find('bsa').text)
                    solv_en = float(residue.find('solv_en').text)
                    if bsa != 0.0:
                        sa_writer.writerow([chain_id, res_name, seq_num, f"{asa:.2f}", f"{bsa:.2f}", f"{solv_en:.2f}"])
                        if chain_id == 'A' and seq_num >= residue_counter:
                            total_bsa += int((bsa / asa) * 100 / 10)
                
                for bond in elem.findall('.//salt-bridges/bond'):
                    sb_writer.writerow([
                        bond.find('chain-1').text,
                        bond.find('res-1').text,
                        bond.find('seqnum-1').text,
                        bond.find('atname-1').text,
                        bond.find('chain-2').text,
                        bond.find('res-2').text,
                        bond.find('seqnum-2').text,
                        bond.find('atname-2').text,
                        "{:.2f}".format(float(bond.find('dist').text))
                    ])
                    chain1 = bond.find('chain-1').text
                    seqnum1 = int(bond.find('seqnum-1').text)
                    chain2 = bond.find('chain-2').text
                    seqnum2 = int(bond.find('seqnum-2').text)
                    if (chain1 == 'A' and seqnum1 > residue_counter) or (chain2 == 'A' and seqnum2 > residue_counter):
                        total_salt_bridges += 1
                
                for bond in elem.findall('.//h-bonds/bond'):
                    hb_writer.writerow([
                        bond.find('chain-1').text,
                        bond.find('res-1').text,
                        bond.find('seqnum-1').text,
                        bond.find('atname-1').text,
                        bond.find('chain-2').text,
                        bond.find('res-2').text,
                        bond.find('seqnum-2').text,
                        bond.find('atname-2').text,
                        "{:.2f}".format(float(bond.find('dist').text))
                    ])
                    chain1 = bond.find('chain-1').text
                    seqnum1 = int(bond.find('seqnum-1').text)
                    chain2 = bond.find('chain-2').text
                    seqnum2 = int(bond.find('seqnum-2').text)
                    if (chain1 == 'A' and seqnum1 > residue_counter) or (chain2 == 'A' and seqnum2 > residue_counter):
                        total_hbonds += 1
            root.clear()
    
    return (base_filename, total_bsa, total_salt_bridges, total_hbonds)

def process_files_in_batches(files, batch_size):
    with open(contacts_file, mode='w', newline='') as file:
        writer = csv.writer(file)
        writer.writerow(header)
        
        for i in range(0, len(files), batch_size):
            batch = files[i:i+batch_size]
            with ThreadPoolExecutor(max_workers=max_workers) as executor:
                futures = {executor.submit(process_file, file): file for file in batch}
                for future in as_completed(futures):
                    result = future.result()
                    if result:
                        writer.writerow(result)

files_to_process = [file for file in os.listdir(INPUT_DIR) if file.startswith("rfdiff_design_") and file.endswith("_af2pred.pdb")]
process_files_in_batches(files_to_process, batch_size)

# Merge contacts and binders data
binders_file = os.path.join(INPUT_DIR, "binders_list.txt")
output_file = os.path.join(OUTPUT_DIR, "final_binders_list.csv")

contacts_exists = os.path.isfile(contacts_file)
binders_exists = os.path.isfile(binders_file)

if not contacts_exists or not binders_exists:
    if not contacts_exists:
        print(f"Error: {contacts_file} not found!")
    if not binders_exists:
        print(f"Error: {binders_file} not found!")
    raise FileNotFoundError("Error: Required files not found!")

with open(contacts_file, mode='r') as contacts, open(binders_file, mode='r') as binders, open(output_file, mode='w', newline='') as output:
    contacts_reader = csv.reader(contacts)
    binders_reader = csv.reader(binders, delimiter='\t')
    output_writer = csv.writer(output)
    
    contacts_data = {rows[0]: rows[1:] for rows in contacts_reader}
    
    binders_header = next(binders_reader)
    output_writer.writerow(binders_header + ["bsa_score", "salt_bridges", "h_bonds"])
    
    for row in binders_reader:
        desc = row[4]
        output_writer.writerow(row + contacts_data.get(desc, ["NA", "NA", "NA"]))

print(f"Merged file has been created: {output_file}")

# Move generated files to tables directory, excluding final_binders_list.csv
tables_dir = os.path.join(OUTPUT_DIR, "tables")
os.makedirs(tables_dir, exist_ok=True)
for file in os.listdir(OUTPUT_DIR):
    file_path = os.path.join(OUTPUT_DIR, file)
    if file.endswith((".csv", ".xml")) and file_path != output_file:
        os.rename(file_path, os.path.join(tables_dir, file))


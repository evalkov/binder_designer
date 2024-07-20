#!/bin/bash
#SBATCH --job-name=contacts
#SBATCH --nodes=1
#SBATCH --cpus-per-task=1
#SBATCH --mem-per-cpu=1G
#SBATCH --time=00:59:00
#SBATCH --mail-type=ALL
#SBATCH --mail-user=$USER

module load ccp4/8.0.010

# This must be set for pisa to work and for pisa.cfg to be read
export CCP4=/mnt/nasapps/production/ccp4/8.0.010

PROC_DIR=/path/
cd $PROC_DIR


# Select residue as the cut-off for final stats, assumes binder is chain A
residue_counter=44
# Creates a table with final stats using the specified cut-off
echo "binder,bsa_score,salt_bridges,h_bonds" > contacts.csv

# Loop through each pdb file and run the analysis commands
for file in rfdiff_design_*_dldesign_*_af2pred.pdb; do
    echo "Processing $file"
    # Extract the base filename without the extension
    base_filename="${file%.pdb}"
    # Run PISA analysis and create an XML output file
    pisa "$file" -analyse "$file" > /dev/null 2>&1
    pisa "$file" -xml interfaces > "${base_filename}.xml"
    pisa "$file" -erase > /dev/null 2>&1
    
    # Print the base filename to contacts.csv
    printf "%s" "$base_filename" >> contacts.csv

    # Process the XML file immediately after creation
    xml="${base_filename}.xml"
    
    # Parse XML and extract data into a CSV file
    awk -F, '
    BEGIN {
        FS = "[<>]"
        OFS = ","
        print "Chain,Residue Name,Sequence Number,ASA,BSA,Solvent Energy"
    }
    /<molecule>/ {
        chain_id = ""
    }
    /<chain_id>/ {
        chain_id = $3
    }
    /<residue>/ {
        res_name = ""
        seq_num = ""
        asa = 0
        bsa = 0
        solv_en = 0
    }
    /<name>/ {
        res_name = $3
    }
    /<seq_num>/ {
        seq_num = $3
    }
    /<asa>/ {
        asa = sprintf("%.2f", $3)
    }
    /<bsa>/ {
        bsa = sprintf("%.2f", $3)
    }
    /<solv_en>/ {
        solv_en = sprintf("%.2f", $3)
    }
    /<\/residue>/ {
        if (bsa != "0.00") {
            print chain_id, res_name, seq_num, asa, bsa, solv_en
        }
    }' "$xml" > "${base_filename}_sa.csv"

# Calculate buried area percentage and create a new CSV with this data
awk -F, '
NR==1 {
    print $0",Buried Area Percentage (%),Buried Area Visual"
    next
}
NR>1 {
    asa = $4
    bsa = $5
    percentage = (bsa / asa) * 100
    bars = int(percentage / 10)
    bar_str = ""
    for (i = 0; i < bars; i++) bar_str = bar_str "|"
    printf "%s,%.2f,%s\n", $0, percentage, bar_str
}' "${base_filename}_sa.csv" > "${base_filename}_sa_with_buried_area_percentage.csv"

# Calculate the total BSA score for residues in chain A with sequence number >= residue_counter
awk -F, -v residue_counter="$residue_counter" '
NR > 1 {
    chain = $1
    seqnum = $3
    bars = length($NF)  # $NF refers to the last field, which is the bar_str
    if (chain == "A" && seqnum >= residue_counter) {
        total_bars += bars
    }
}
END {
    print "Total BSA score is:", total_bars
    printf ",%s", total_bars >> "contacts.csv"
}' "${base_filename}_sa_with_buried_area_percentage.csv"


# Parse XML and extract salt bridge data into a CSV file
awk -F, '
BEGIN {
    FS = "[<>]"
    OFS = ","
    print "Chain 1,Residue 1,Sequence Number 1,Atom Name 1,Chain 2,Residue 2,Sequence Number 2,Atom Name 2,Distance"
}

/<molecule>/ {
    chain_id_1 = ""
    chain_id_2 = ""
}

/<chain_id>/ {
    if (chain_id_1 == "") {
        chain_id_1 = $3
    } else {
        chain_id_2 = $3
    }
}

/<salt-bridges>/ {
    capturing_salt_bridges = 1
}

/<\/salt-bridges>/ {
    capturing_salt_bridges = 0
}

capturing_salt_bridges && /<bond>/ {
    bond_chain_1 = ""
    bond_res_1 = ""
    bond_seqnum_1 = ""
    bond_atname_1 = ""
    bond_chain_2 = ""
    bond_res_2 = ""
    bond_seqnum_2 = ""
    bond_atname_2 = ""
    bond_dist = ""
}

capturing_salt_bridges && /<chain-1>/ {
    bond_chain_1 = $3
}

capturing_salt_bridges && /<res-1>/ {
    bond_res_1 = $3
}

capturing_salt_bridges && /<seqnum-1>/ {
    bond_seqnum_1 = $3
}

capturing_salt_bridges && /<atname-1>/ {
    bond_atname_1 = $3
}

capturing_salt_bridges && /<chain-2>/ {
    bond_chain_2 = $3
}

capturing_salt_bridges && /<res-2>/ {
    bond_res_2 = $3
}

capturing_salt_bridges && /<seqnum-2>/ {
    bond_seqnum_2 = $3
}

capturing_salt_bridges && /<atname-2>/ {
    bond_atname_2 = $3
}

capturing_salt_bridges && /<dist>/ {
    bond_dist = sprintf("%.2f", $3)
}

capturing_salt_bridges && /<\/bond>/ {
    print bond_chain_1, bond_res_1, bond_seqnum_1, bond_atname_1, bond_chain_2, bond_res_2, bond_seqnum_2, bond_atname_2, bond_dist
}' "${base_filename}.xml" > "${base_filename}_salt_bridges.csv"

# Calculate the total number of salt bridges involving residues in chain A with sequence number > residue_counter
awk -F, -v residue_counter="$residue_counter" '
BEGIN {
    total_salt_bridges = 0
}
NR > 1 {
    chain1 = $1
    seqnum1 = $3
    chain2 = $5
    seqnum2 = $7
    if ((chain1 == "A" && seqnum1 > residue_counter) || (chain2 == "A" && seqnum2 > residue_counter)) {
        total_salt_bridges += 1
    }
}
END {
    print "Total number of salt bridges is:", total_salt_bridges
    printf ",%s", total_salt_bridges >> "contacts.csv"
}' "${base_filename}_salt_bridges.csv"

# Parse XML and extract hydrogen bond data into a CSV file
awk -F, '
BEGIN {
    FS = "[<>]"
    OFS = ","
    print "Chain 1,Residue 1,Sequence Number 1,Atom Name 1,Chain 2,Residue 2,Sequence Number 2,Atom Name 2,Distance"
}

/<molecule>/ {
    chain_id_1 = ""
    chain_id_2 = ""
}

/<chain_id>/ {
    if (chain_id_1 == "") {
        chain_id_1 = $3
    } else {
        chain_id_2 = $3
    }
}

/<h-bonds>/ {
    capturing_h_bonds = 1
}

/<\/h-bonds>/ {
    capturing_h_bonds = 0
}

capturing_h_bonds && /<bond>/ {
    bond_chain_1 = ""
    bond_res_1 = ""
    bond_seqnum_1 = ""
    bond_atname_1 = ""
    bond_chain_2 = ""
    bond_res_2 = ""
    bond_seqnum_2 = ""
    bond_atname_2 = ""
    bond_dist = ""
}

capturing_h_bonds && /<chain-1>/ {
    bond_chain_1 = $3
}

capturing_h_bonds && /<res-1>/ {
    bond_res_1 = $3
}

capturing_h_bonds && /<seqnum-1>/ {
    bond_seqnum_1 = $3
}

capturing_h_bonds && /<atname-1>/ {
    bond_atname_1 = $3
}

capturing_h_bonds && /<chain-2>/ {
    bond_chain_2 = $3
}

capturing_h_bonds && /<res-2>/ {
    bond_res_2 = $3
}

capturing_h_bonds && /<seqnum-2>/ {
    bond_seqnum_2 = $3
}

capturing_h_bonds && /<atname-2>/ {
    bond_atname_2 = $3
}

capturing_h_bonds && /<dist>/ {
    bond_dist = sprintf("%.2f", $3)
}

capturing_h_bonds && /<\/bond>/ {
    print bond_chain_1, bond_res_1, bond_seqnum_1, bond_atname_1, bond_chain_2, bond_res_2, bond_seqnum_2, bond_atname_2, bond_dist
}' "${base_filename}.xml" > "${base_filename}_hbonds.csv"

# Calculate the total number of hydrogen bonds involving residues in chain A with sequence number > residue_counter
awk -F, -v residue_counter="$residue_counter" '
BEGIN {
    total_hbonds = 0
}
NR > 1 {
    chain1 = $1
    seqnum1 = $3
    chain2 = $5
    seqnum2 = $7
    if ((chain1 == "A" && seqnum1 > residue_counter) || (chain2 == "A" && seqnum2 > residue_counter)) {
        total_hbonds += 1
    }
}
END {
    printf ",%s\n", total_hbonds >> "contacts.csv"
    print "Total number of H bonds is:", total_hbonds
}' "${base_filename}_hbonds.csv"

done

mkdir -p tables
mv *_sa.csv tables/
mv *_sa_with_buried_area_percentage.csv tables/
mv *_salt_bridges.csv tables/
mv *_hbonds.csv tables/
mv *.xml tables/

# Define file paths
contacts_file="contacts.csv"
binders_file="binders_list.txt"
output_file="final_binders_list.csv"

# Check if required files exist
if [[ ! -f "$contacts_file" ]] || [[ ! -f "$binders_file" ]]; then
  echo "Error: Required files not found!"
  exit 1
fi

# Extract headers from both files
binders_header=$(head -n 1 "$binders_file")
contacts_header=$(head -n 1 "$contacts_file")

# Create the output file and add the merged header
echo "${binders_header},bsa_score,salt_bridges,h_bonds" | tr '\t' ',' > "$output_file"

# Read contacts.csv into a temporary file
awk -F, 'NR > 1 {print $1 "," $2 "," $3 "," $4}' "$contacts_file" > contacts_temp.txt

# Read binders_list.txt and join with contacts data
awk -F'\t' '
  BEGIN {
    OFS = ","
    while ((getline line < "contacts_temp.txt") > 0) {
      split(line, fields, ",")
      contacts[fields[1]] = fields[2] "," fields[3] "," fields[4]
    }
  }
  NR > 1 {
    desc = $5
    if (desc in contacts) {
      print $0 OFS contacts[desc]
    } else {
      print $0 OFS "NA,NA,NA"
    }
  }
' "$binders_file" | tr '\t' ',' >> "$output_file"

# Clean up the temporary file
rm contacts_temp.txt

echo "Merged file has been created: $output_file"

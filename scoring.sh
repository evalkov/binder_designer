#!/bin/bash

# Check if RFDIR is set
if [ -z "$RFDIR" ]; then
  echo "RFDIR is not set. Exiting."
  exit 1
fi

cd $RFDIR

AF2_PAE_INTERACT="10"

# Check if the file out.sc exists in the current directory or in the logs directory
if [[ -e "out.sc" ]]; then
    input_file="out.sc"
    echo "File out.sc exists in the current directory."
elif [[ -e "logs/out.sc" ]]; then
    ln -s logs/out.sc out.sc
    input_file="out.sc"
    echo "File out.sc does not exist in both the current directory and logs."
    # Handle the case where the file doesn't exist, e.g., exit or set a default
fi
# Process the out.sc file
{
    # Extract and process the header: remove the first column
    head -n 1 "$input_file" | awk '{$1=""; print $0}' OFS="\t" | sed 's/^\t//'  | awk 'BEGIN {OFS="\t"} {print $1, $2, $3, $5, $10}'

    # Process the data: convert to tab-delimited, skip the first line, remove the first column,
    # filter for pae_interaction < 10, and sort by pae_interaction
    awk '{$1=$1; print}' OFS="\t" "$input_file" | cut -f2- | awk -F'\t' '$3 < '$AF2_PAE_INTERACT'' | sort -k3,3n | awk 'BEGIN {OFS="\t"} {print $1, $2, $3, $5, $10}'
} > pdb_list.txt


# Check if binders identified, and if not, then send an email and exit
if ! grep -q "_af2pred" pdb_list.txt; then
    echo "No designed binders had acceptable metrics (pae_interaction <10)." | mutt -s "RFdiffusion complete." -e 'my_hdr From:RFdiffusion (RFdiffusion)' -- $USER@nih.gov
    exit 1
fi

sed '1s/$/\tbinder_seq\tbinder_len/' pdb_list.txt | head -n 1 > binders_list.txt

# Read each line from pdb_list.txt, skipping the first header line
tail -n +2 pdb_list.txt | while read -r line; do
    # Extract the filename from the last column
    filename=$(echo "$line" | awk '{print $NF}').pdb
    list_chain() {
                awk '/^ATOM/ && $3 == "CA" {print $5}' | uniq
        }
        extract_seq_chain() {
                awk -v ch=$1 '/^ATOM/ && $3 == "CA" && $5 == ch {print $4}'
        }
        remove_newline() {
                tr '\n' ' '
        }
        convert_aa() {
                sed 's/ALA/A/g;s/CYS/C/g;s/ASP/D/g;s/GLU/E/g;s/PHE/F/g;s/GLY/G/g;s/HIS/H/g;s/ILE/I/g;s/LYS/K/g;s/LEU/L/g;s/MET/M/g;s/ASN/N/g;s/PRO/P/g;s/GLN/Q/g;s/ARG/R/g;s/SER/S/g;s/THR/T/g;s/VAL/V/g;s/TRP/W/g;s/TYR/Y/g'
        }
        remove_space() {
                sed 's/ //g'
        }
    # Check if the file exists
    if [[ -e "$filename" ]]; then
        sequence=$(cat "$filename" | extract_seq_chain 'A' | remove_newline | convert_aa | remove_space)
    elif [[ -e "af2_models/$filename" ]]; then
        sequence=$(cat "af2_models/$filename" | extract_seq_chain 'A' | remove_newline | convert_aa | remove_space)
    fi
    if [[ -n $sequence ]]; then
                size=$(echo $sequence | wc -c)
                size=$((size-1))
                echo -e "$line\t$sequence\t$size" | tee -a binders_list.txt
    fi
done

mkdir -p top_scoring_solutions

# Process lines starting from line 2
awk 'NR > 1 {print $5 ".pdb"}' binders_list.txt | while read -r file; do
    # Check if the file exists
    if [[ -e $file ]]; then
        cp "$file" top_scoring_solutions/
    fi
done

mkdir -p top20_scoring_solutions
cp binders_list.txt top20_scoring_solutions/
cp binders_list.txt top_scoring_solutions/

echo "\
set bgcolor white" > binders.cxc

# Process lines starting from line 2 and copy first 20 files
awk 'NR > 1 && NR <= 21 {print $5 ".pdb"}' binders_list.txt | while read -r file; do
    # Check if the file exists
    if [[ -e $file ]]; then
        echo "open $file" >> binders.cxc
        cp "$file" top20_scoring_solutions/
    elif [[ -e af2_models/$file ]]; then
        echo "open $file" >> binders.cxc
        cp af2_models/"$file" top20_scoring_solutions/
    fi
done

echo "\
cartoon style protein modeh tube rad 2 sides 24
cartoon style width 2 thick 0.2
rainbow chain palette RdYlBu-5
lighting simple shadows false intensity 0.5
view all
hide atoms
show cartoons
hide all models
show #1 models
matchmaker all to #1/B pairing bs
" >> binders.cxc

# Initialize a counter
counter=0

# Loop through each line of the input file, skipping the header line
awk 'NR > 1 && NR <= 21 {print $5 ".pdb"}' binders_list.txt | while read -r file; do
    # Increment the counter
    ((counter++))

    # Check if the counter has reached 20, if so, exit the loop
    if [[ $counter -eq 21 ]]; then
        break
    fi

    # Check if the file exists
    if [[ -e $file ]]; then
        echo "\
interfaces select #$counter/B contacting #$counter/A bothSides true
contacts #$counter/A restrict #$counter/B intraMol false
show sel atoms
select clear
" >> binders.cxc
    elif [[ -e af2_models/$file ]]; then
        echo "\
interfaces select #$counter/B contacting #$counter/A bothSides true
contacts #$counter/A restrict #$counter/B intraMol false
show sel atoms
select clear
" >> binders.cxc
    fi
done

echo "\
delete H
color byhetero
" >> binders.cxc

cp binders.cxc top20_scoring_solutions/
cp binders.cxc top_scoring_solutions/

# Count lines excluding the first line in each file
hit_binders=$(tail -n +2 binders_list.txt | wc -l)

tar -cvjf top20_scoring_solutions.tar.bz2 top20_scoring_solutions/
tar -cvjf top_scoring_solutions.tar.bz2 top_scoring_solutions/

# Check files size before attachign to email
maxsize=$((7*1024*1024)) # 7 MB in bytes
totalsize=$(stat -c%s top20_scoring_solutions.tar.bz2)

# workaround with libcrypto not accessing the correct openssl
export LD_LIBRARY_PATH=/lib64:\$LD_LIBRARY_PATH

if [ $totalsize -lt $maxsize ]; then
        cat rfdiff.inp | mutt -e 'set content_type=text/html' -s "RFdiffusion is complete with $hit_binders hits." -e 'my_hdr From:RFdiffusion (RFdiffusion)' -a top20_scoring_solutions.tar.bz2 -- $USER@nih.gov
else
        cat rfdiff.inp | mutt -e 'set content_type=text/html' -s "RFdiffusion is complete with $hit_binders hits." -e 'my_hdr From:RFdiffusion (RFdiffusion)' -- $USER@nih.gov
fi

#!/bin/bash
#SBATCH --job-name=scoring
#SBATCH --cpus-per-task=48
#SBATCH --mem-per-cpu=1G
#SBATCH --time=00:30:00
#SBATCH --mail-type=FAIL,INVALID_DEPEND,TIME_LIMIT

module load python
module load ccp4/8.0.010

# Define the script directory
SOFT=/home/valkove2/soft/binder_designer

# This must be set for pisa to work and for pisa.cfg to be read
export CCP4=/mnt/nasapps/production/ccp4/8.0.010

AF2_PAE_INTERACT="10"

if [[ -e "out_1.sc" ]]; then
    cp out_1.sc out.sc
    for file in $(ls out_*.sc); do
        if [ "$file" != "out_1.sc" ]; then
            tail -n +2 $file >> out.sc
        fi
    done
else
    echo "out_1.sc doesn't exist"
    exit 1
fi

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

# Check if a file exists and exit if it doesn't
check_file_exists() {
    if [ ! -f "$1" ]; then
        echo "File $1 does not exist. Exiting."
        exit 1
    fi
}

# Execute a command and exit if it fails
execute_command() {
    eval $1
    if [ $? -ne 0 ]; then
        echo "Command failed: $1. Exiting."
        exit 1
    fi
}

# Run the Python scripts with the appropriate checks
execute_command "`which python` $SOFT/compute.py"
execute_command "`which python` $SOFT/pisa.py"
execute_command "`which python` $SOFT/binder_list_generation.py"
execute_command "`which python` $SOFT/scaler.py"

# Check if input files exist before running heatmap generation scripts
check_file_exists "top_50_binders_weighted_minmax_scaler.csv"
execute_command "`which python` $SOFT/heatmap_generation.py top_50_binders_weighted_minmax_scaler.csv"

check_file_exists "top_50_binders_weighted_standard_scaler.csv"
execute_command "`which python` $SOFT/heatmap_generation.py top_50_binders_weighted_standard_scaler.csv"

check_file_exists "binders_weighted_composite_scores_minmax_scaler.csv"
execute_command "`which python` $SOFT/plots_generation.py binders_weighted_composite_scores_minmax_scaler.csv"

check_file_exists "binders_weighted_composite_scores_standard_scaler.csv"
execute_command "`which python` $SOFT/plots_generation.py binders_weighted_composite_scores_standard_scaler.csv"

# Check if both input files exist before running the compare_csv script
check_file_exists "top_50_binders_weighted_standard_scaler.csv"
check_file_exists "top_50_binders_weighted_minmax_scaler.csv"
execute_command "`which python` $SOFT/compare_csv.py top_50_binders_weighted_standard_scaler.csv top_50_binders_weighted_minmax_scaler.csv"
execute_command "`which python` $SOFT/top_binder_scatter_plot_generation.py"
execute_command "`which python` $SOFT/chimerax_script_generation.py"

# Count lines excluding the first line in each file
hit_binders=$(tail -n +2 binders_list.txt | wc -l)

# Check files size before attachign to email
maxsize=$((7*1024*1024)) # 7 MB in bytes
totalsize=$(stat -c%s predictions.tar.bz2)

# workaround with libcrypto not accessing the correct openssl
export LD_LIBRARY_PATH=/lib64:\$LD_LIBRARY_PATH

if [ $totalsize -lt $maxsize ]; then
        cat run_parameters.html | mutt -e 'set content_type=text/html' -s "RFdiffusion is complete with $hit_binders hits." -e 'my_hdr From:RFdiffusion (RFdiffusion)' -a predictions.tar.bz2 -- $USER@nih.gov
else
        cat run_parameters.html | mutt -e 'set content_type=text/html' -s "RFdiffusion is complete with $hit_binders hits." -e 'my_hdr From:RFdiffusion (RFdiffusion)' -- $USER@nih.gov
fi

# Tidy up
# Create directories
mkdir -p scoring_tables plots logs af2pred_pdbs silent_files af2_scores temp_files processing

# Move files to respective directories
mv *.csv scoring_tables/
mv *.log logs/
mv *.eps plots/
mv *_af2pred.pdb af2pred_pdbs/
mv *.silent silent_files/
mv *.sc af2_scores/
mv *.point *.idx x?? temp_files/
mv rfdiff_chunk_* silent_files pisa_xml_files temp_files af2_scores pdb trb processing/

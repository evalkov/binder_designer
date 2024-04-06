#!/bin/bash
#SBATCH --job-name=rfdiff
#SBATCH --partition=gpu
#SBATCH --gres=gpu:v100:2
#SBATCH --nodes=1
#SBATCH --cpus-per-task=16
#SBATCH --mem-per-cpu=1G
#SBATCH --time=4-00:00:00
#SBATCH --mail-type=ALL
#SBATCH --mail-user=valkove2

model="/home/valkove2/pdbs/dnd1.pdb"
binder_length="70-100"
hotspot_residue="57"
num_binders="1000"

RFDIFFUSION_DIR="/home/valkove2/soft/RFdiffusion"
OUT_DIR="/mnt/beegfs/valkov/rfdiffusion"
DL_BINDER_DESIGN_DIR="/home/valkove2/soft/dl_binder_design"
SILENT_TOOLS="/home/valkove2/soft/silent_tools"

model_id=`basename "$model" .pdb`
hotspot=`awk '/^ATOM/ && $6 == '$hotspot_residue' {residue_name = $4} END {print residue_name}' $model`
template_chain=`awk '/^ATOM/ {chain_id = substr($0, 22, 1)} END {print chain_id}' $model`
template_length=`awk '/^ATOM/ {residue = $6} END {print residue}' $model`
target=$(sed 's|.*/||; s|\.pdb$||' <<< "$model")
timestamp=$(date +"%Y%m%d_%H%M%S")
project=$(sed 's|.*/||; s|\.pdb$||' <<< "$model")_"$timestamp"
mkdir $OUT_DIR/"$project"
cp $model $OUT_DIR/$project

module purge
module load rfdiffusion/1.1.0

$RFDIFFUSION_DIR/scripts/run_inference.py \
	inference.output_prefix=$OUT_DIR/$project/rfdiff \
	inference.input_pdb=$model \
	'contigmap.contigs=['$template_chain''1'-'$template_length'/0 '$binder_length']' \
	'ppi.hotspot_res=['$template_chain''$hotspot_residue']' \
	inference.num_designs=$num_binders \
	inference.ckpt_override_path=$RFDIFFUSION_DIR/models/Complex_beta_ckpt.pt

unset PYTHONPATH
module purge
export PATH=${SILENT_TOOLS}:$PATH
module load dl_binder_design/proteinmpnn_binder_design

cd $OUT_DIR/$project

silentfrompdbs rfdiff_*.pdb > data.silent

$DL_BINDER_DESIGN_DIR/mpnn_fr/dl_interface_design.py -silent data.silent -outsilent mpnn.silent

module purge
module load dl_binder_design/af2_binder_design
module load tensorRT/8.6.1-cuda12 cuda/12.0

$DL_BINDER_DESIGN_DIR/af2_initial_guess/predict.py -silent mpnn.silent -outsilent af2.silent

silentextract af2.silent

# Define the input file name
input_file="out.sc"

# Process the file
{
    # Extract and process the header: remove the first column
    head -n 1 "$input_file" | awk '{$1=""; print $0}' OFS="\t" | sed 's/^\t//'

    # Process the data: convert to tab-delimited, skip the first line, remove the first column,
    # filter for pae_interaction < 10, and sort by pae_interaction
    awk '{$1=$1; print}' OFS="\t" "$input_file" | tail -n +2 | cut -f2- | awk -F'\t' '$3 < 10' | sort -k3,3n
} > filtered_sorted_output.txt

# Display the filtered and sorted data
echo -e "\
Top scoring solutions based on pae_interaction < 10:\n"
cat filtered_sorted_output.txt

mkdir top_scoring_solutions

awk 'NR > 1 {print $NF ".pdb"}' filtered_sorted_output.txt | while IFS= read -r file; do
        cp "$file" top_scoring_solutions/
done

for file in top_scoring_solutions/*.pdb; do
	echo "open $(basename "$file")" >> binders.cxc
done

echo "\
set bgcolor white
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

mv binders.cxc top_scoring_solutions/

# Define the output file
output_file="top_scoring_solutions.html"

# Start the HTML table and write to the output file
echo '<table style="border-collapse: collapse; width: 100%;">' > "$output_file"

# Initialize the variable to check for the first line
first_line=true

# Read data line by line from the file
while read -r line; do
    # Check if it's the first line (header)
    if [[ $first_line == true ]]; then
        echo '  <tr>' >> "$output_file"
        for word in $line; do
            echo "    <th style=\"border: 1px solid black; padding: 8px;\">$word</th>" >> "$output_file"
        done
        echo '  </tr>' >> "$output_file"
        first_line=false
    else
        echo '  <tr>' >> "$output_file"
        for word in $line; do
            echo "    <td style=\"border: 1px solid black; padding: 8px;\">$word</td>" >> "$output_file"
        done
        echo '  </tr>' >> "$output_file"
    fi
done < filtered_sorted_output.txt

# Close the table and write to the output file
echo '</table>' >> "$output_file"

# Count lines excluding the first line in each file
total_binders=$(tail -n +2 out.sc | wc -l)
hit_binders=$(tail -n +2 filtered_sorted_output.txt | wc -l)

percent_hits=$(echo "scale=2; $hit_binders / $total_binders * 100" | bc)

cp top_scoring_solutions.html top_scoring_solutions/

tar -cvjf top_scoring_solutions.tar.bz2 top_scoring_solutions/

# Check files size before attachign to email
maxsize=$((7*1024*1024)) # 7 MB in bytes
size=$(stat -c%s top_scoring_solutions.tar.bz2)

# workaround with libcrypto not accessing the correct openssl
export LD_LIBRARY_PATH=/lib64:\$LD_LIBRARY_PATH

if [ $size -lt $maxsize ]; then
	cat top_scoring_solutions.html | mutt -e 'set content_type=text/html' -s "Binder design for $target is complete with $percent_hits% hits." -e 'my_hdr From:RFdiffusion (RFdiffusion)' -a top_scoring_solutions.tar.bz2 -- $USER@nih.gov
else
	cat top_scoring_solutions.html | mutt -e 'set content_type=text/html' -s "Binder design for $target is complete with $percent_hits% hits." -e 'my_hdr From:RFdiffusion (RFdiffusion)' -- $USER@nih.gov
fi

# Clean up
rm filtered_sorted_output.txt
mkdir rfdiff_models
mkdir af2_models
mkdir silent_files
mkdir logs
mv rfdiff_*.pdb rfdiff_models/
mv rfdiff_models/rfdiff_*_af2pred.pdb af2_models/
mv rfdiff_*.trb rfdiff_models/
mv *.silent* silent_files/
mv out.sc logs/
mv check.point logs/


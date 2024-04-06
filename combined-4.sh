#!/bin/bash
#SBATCH --job-name=BEENDING
#SBATCH --partition=gpu
#SBATCH --gres=gpu:v100:2
#SBATCH --nodes=1
#SBATCH --cpus-per-task=16
#SBATCH --mem-per-cpu=1G
#SBATCH --time=4-00:00:00
#SBATCH --mail-type=ALL
#SBATCH --mail-user=valkove2

model="/home/valkove2/pdbs/not9.pdb"
binder_length="20-30"
hotspot_residue="166"
num_binders="5"
betamodel="no"

RFDIFFUSION_DIR="/home/valkove2/soft/RFdiffusion"
OUT_DIR="/mnt/beegfs/valkov/rfdiffusion"
DL_BINDER_DESIGN_DIR="/home/valkove2/soft/dl_binder_design"
SILENT_TOOLS="/home/valkove2/soft/silent_tools"

model_id=`basename $model`
hotspot=`awk '/^ATOM/ && $6 == '$hotspot_residue' {residue_name = $4} END {print residue_name}' $model`
template_chain=`awk '/^ATOM/ {chain_id = substr($0, 22, 1)} END {print chain_id}' $model`
template_length=`awk '/^ATOM/ {residue = $6} END {print residue}' $model`
timestamp=$(date +"%Y%m%d_%H%M%S")
project=$(sed 's|.*/||; s|\.pdb$||' <<< "$model")_"$timestamp"

mkdir $OUT_DIR/"$project"
cp $model $OUT_DIR/$project
cd $OUT_DIR/$project

echo "\
BINDER DESIGN PIPELINE 
Written by Eugene Valkov, NCI/NIH.
" | tee -a $OUT_DIR/$project/settings.txt

if [ "$betamodel" == "yes" ]; then
        echo -e "\
RFdiffusion v `awk -F"[',]" '/version=/{print $2}' $RFDIFFUSION_DIR/setup.py` with Complex_beta_ckpt.pt model weight used to design" | tee -a $OUT_DIR/$project/settings.txt
else
	echo -e "\
RFdiffusion v `awk -F"[',]" '/version=/{print $2}' $RFDIFFUSION_DIR/setup.py` with default model weight used to design" | tee -a $OUT_DIR/$project/settings.txt
fi

echo "\
$num_binders binders of $binder_length residues targeting $hotspot-$hotspot_residue of $model_id, chain $template_chain ($template_length residues)
Sequence assigned with ProteinMPNN and validated with AlphaFold via dl_binder_design
Computed on `hostname`
All files at $OUT_DIR/$project
Started on on `date "+%B %d, %Y"` at `date "+%I:%M %p"`" | tee -a $OUT_DIR/$project/settings.txt

module purge
module load rfdiffusion/1.1.0

# Construct the command with mandatory arguments
cmd="$RFDIFFUSION_DIR/scripts/run_inference.py \
        inference.output_prefix=$OUT_DIR/$project/rfdiff \
        inference.input_pdb=$model \
        'contigmap.contigs=['$template_chain''1'-'$template_length'/0 '$binder_length']' \
        'ppi.hotspot_res=['$template_chain''$hotspot_residue']' \
        inference.num_designs=$num_binders"

# Append the optional argument if the condition is met
if [ "$betamodel" == "yes" ]; then
    cmd+=" inference.ckpt_override_path=$RFDIFFUSION_DIR/models/Complex_beta_ckpt.pt"
fi

# Execute the command
eval $cmd | tee -a $OUT_DIR/$project/rfdiffusion.log

unset PYTHONPATH
module purge
export PATH=${SILENT_TOOLS}:$PATH
module load dl_binder_design/proteinmpnn_binder_design

#cd $OUT_DIR/$project

silentfrompdbs rfdiff_*.pdb > data.silent

$DL_BINDER_DESIGN_DIR/mpnn_fr/dl_interface_design.py -silent data.silent -outsilent mpnn.silent | tee -a $OUT_DIR/$project/proteinmpnn.log

module purge
module load dl_binder_design/af2_binder_design
module load tensorRT/8.6.1-cuda12 cuda/12.0

$DL_BINDER_DESIGN_DIR/af2_initial_guess/predict.py -silent mpnn.silent -outsilent af2.silent | tee -a $OUT_DIR/$project/alphafold.log

silentextract af2.silent

input_file=out.sc

# Process the file
{
    # Extract and process the header: remove the first column
    head -n 1 "$input_file" | awk '{$1=""; print $0}' OFS="\t" | sed 's/^\t//'  | awk 'BEGIN {OFS="\t"} {print $1, $2, $3, $5, $10}'

    # Process the data: convert to tab-delimited, skip the first line, remove the first column,
    # filter for pae_interaction < 10, and sort by pae_interaction
    awk '{$1=$1; print}' OFS="\t" "$input_file" | cut -f2- | awk -F'\t' '$3 < 10' | sort -k3,3n | awk 'BEGIN {OFS="\t"} {print $1, $2, $3, $5, $10}'
} > pdb_list.txt


# Check if binders identified, and if not, then send an email and exit
if ! grep -q "_dldesign_0_cycle1_af2pred" pdb_list.txt; then
    echo "\
Completed on `date "+%B %d, %Y"` at `date "+%I:%M %p"`
" | tee -a $OUT_DIR/$project/settings.txt
    echo -e "\
No designed binders had acceptable metrics (pae_interaction <10). 
Try to increase the number of initial designs, change the binder chain length, or choose a different hotspot." | tee -a settings.txt
    cat settings.txt | mutt -s "Binder design for $model_id found no hits." -e 'my_hdr From:RFdiffusion (RFdiffusion)' -- $USER@nih.gov
    exit 1
fi

sed '1s/$/\tbinder_seq\tbinder_len/' pdb_list.txt | head -n 1 > binders_list.txt

# Read each line from pdb_list.txt, skipping the first header line
tail -n +2 pdb_list.txt | while read -r line; do
    # Extract the filename from the last column
    filename=$(echo "$line" | awk '{print $NF}').pdb
    # Check if the file exists
    if [[ -e "$filename" ]]; then
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
	sequence=$(cat "$filename" | extract_seq_chain 'A' | remove_newline | convert_aa | remove_space)
	if [[ -n $sequence ]]; then
        	size=$(echo $sequence | wc -c)
        	size=$((size-1))
		echo -e "$line\t$sequence\t$size" | tee -a binders_list.txt
    	fi
    fi
done

echo "\
Completed on `date "+%B %d, %Y"` at `date "+%I:%M %p"`

Top-scoring binder designs are shown on the right.
" | tee -a $OUT_DIR/$project/settings.txt

mkdir top_scoring_solutions

# Process lines starting from line 2
awk 'NR > 1 {print $5 ".pdb"}' binders_list.txt | while read -r file; do
    # Check if the file exists
    if [[ -e $file ]]; then
        cp "$file" top_scoring_solutions/
    fi
done

mkdir top20_scoring_solutions
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
    fi
done

echo "\
delete H
color byhetero
" >> binders.cxc

cp binders.cxc top20_scoring_solutions/
cp binders.cxc top_scoring_solutions/

# Make a movie of top hits

echo "open $model_id" > movie.cxc

# Loop through each line of the input file, skipping the header line
awk 'NR > 1 && NR <= 11 {print $5 ".pdb"}' binders_list.txt | while read -r file; do
    # Check if the file exists
    if [[ -e $file ]]; then
            echo "open $file" >> movie.cxc
    fi
done

echo "\
set bgcolor white
cartoon style protein modeh tube rad 2 sides 24
cartoon style width 2 thick 0.2
rainbow chain palette RdYlBu-5
lighting flat shadows false intensity 0.3
matchmaker all to #1 pairing bs
select /$template_chain
del H
cartoon hide sel
show #1 cartoons
color #1 cornflower blue
select #1:$hotspot_residue
hide sel atoms
show sel cartoons
show sel atoms
style sel sphere
color sel byelement
select clear
cofr #1:$hotspot_residue
view cofr false
hide all models
show #1 models
movie record size 1500,1500 format png supersample 4 directory $OUT_DIR/$project
" >> movie.cxc

# Initialize a counter
counter=1
# Loop through each line of the input file, skipping the header line
awk 'NR > 1 && NR <= 11 {print $5 ".pdb"}' binders_list.txt | while read -r file; do
    # Increment the counter
    ((counter++))

    # Check if the counter has reached 20, if so, exit the loop
    if [[ $counter -eq 11 ]]; then
        break
    fi

    # Check if the file exists
    if [[ -e $file ]]; then
        echo "\
show #$counter models
turn y 12 360 models #1-$counter
2dlab text '`echo $file | sed 's|.*/||; s|\.pdb$||'`' color black size 10 x .03 y .03
wait 30
stop
hide #$counter models
2dlab delete
" >> movie.cxc
    fi
done

echo "movie stop" >> movie.cxc

module load ChimeraX/1.7

chimerax --offscreen --script movie.cxc --exit

convert \
	-dispose previous \
	-delay 10 \
       	-loop 0 \
       	-dither None \
	-colors 256 \
	-layers Optimize \
	-resize 350x350 \
	-filter Lanczos \
	-coalesce \
	chimovie*.png \
       	animated.gif

cp $model_id top20_scoring_solutions/

# Create an HTML table with top-scorign hits and embedded GIF

# Define the input and output files
input_file="binders_list.txt"
settings_file="settings.txt"
output_file="binders.html"

# Check if the input file exists
if [ ! -f "$input_file" ]; then
    echo "File not found: $input_file"
    exit 1
fi

# Check if the settings file exists
if [ ! -f "$settings_file" ]; then
    echo "File not found: $settings_file"
    exit 1
fi

# Function to include settings
include_settings() {
    echo "<div class=\"settings-title\">"
    # Read each line of settings file and add it to the HTML with a line break.
    while IFS= read -r line; do
        echo "$line<br>"
    done < "$settings_file"
    echo "</div>"
}

# Function to generate HTML table from data in a file
generate_html_table() {
    {
        echo "<!DOCTYPE html>"
        echo "<html>"
        echo "<head>"
        echo "<style>"
        echo "body {"
        echo "  font-family: Arial, sans-serif;"
        echo "}"
        echo ".settings-title {"
        echo "  background-color:  white;"
        echo "  color: black;"
        echo "  font-size: 20px;"
        echo "  font-weight: bold;"
        echo "  text-align: left;"
        echo "  padding: 20px;"
        echo "  margin-top: 20px;"
        echo "}"
        echo "table {"
        echo "  border-collapse: collapse;"
        echo "  width: 100%;"
        echo "  box-shadow: 0 2px 3px #c8d1d3;"
        echo "}"
        echo "th, td {"
        echo "  border: 1px solid #ddd;"
        echo "  padding: 8px;"
        echo "  text-align: left;"
        echo "}"
        echo "th {"
        echo "  background-color: #666;"
        echo "  color: white;"
        echo "}"
        echo "tr:nth-child(even) {"
        echo "  background-color: #f2f2f2;"
        echo "}"
        echo "</style>"
        echo "<style>"
        echo ".header-container {"
        echo "  display: flex;"
        echo "  justify-content: space-between;"
        echo "  align-items: center;"
        echo "}"
        echo ".header-img {"
        echo "  height: 100px;" # Adjust the size of the GIF as needed
        echo "  padding-left: 150px;" # Add some space between the header and the GIF
        echo "}"
        echo "</style>"
        echo "</head>"
        echo "<body>"
        echo "<div class=\"header-container\">"
        echo "  <!-- Existing header content goes here. Make sure this contains your actual header content. -->"
        echo "</div>"
        echo "<div style=\"display: flex; align-items: stretch;\">"

        # Include settings above the table
        include_settings

        echo "  <img alt=\"Animated GIF\" class=\"header-img\" src=\"animated.gif\" style=\"height: 100%; width: auto;\"/>"
        echo "</div>"
        echo "<table>"

        # Read the first line for headers
        read -r headers < "$input_file"
        echo "<tr>"
        for header in $headers; do
            echo "<th>${header}</th>"
        done
        echo "</tr>"

        # Skip the first line and read the rest of the data
        tail -n +2 "$input_file" | while IFS= read -r line; do
            echo "<tr>"
            for value in $line; do
                echo "<td>${value}</td>"
            done
            echo "</tr>"
        done

        echo "</table>"
        echo "</body>"
        echo "</html>"
    } > "$output_file"
}

# Generate the HTML table and output to a file
generate_html_table

echo "HTML table generated in file $output_file"


# Count lines excluding the first line in each file
total_binders=$(tail -n +2 out.sc | wc -l)
hit_binders=$(tail -n +2 binders_list.txt | wc -l)

percent_hits=$(echo "scale=2; $hit_binders / $total_binders * 100" | bc)

tar -cvjf top20_scoring_solutions.tar.bz2 top20_scoring_solutions/
tar -cvjf top_scoring_solutions.tar.bz2 top_scoring_solutions/

# Check files size before attachign to email
maxsize=$((7*1024*1024)) # 7 MB in bytes
size=$(stat -c%s top20_scoring_solutions.tar.bz2)

size1=$(stat -c%s top20_scoring_solutions.tar.bz2)
size2=$(stat -c%s binders.html)
size3=$(stat -c%s animated.gif)
totalsize=$((size1 + size2 + size3))

# workaround with libcrypto not accessing the correct openssl
export LD_LIBRARY_PATH=/lib64:\$LD_LIBRARY_PATH

if [ $totalsize -lt $maxsize ]; then
	cat binders.html | mutt -e 'set content_type=text/html' -s "Binder design for $model_id is complete with $percent_hits% hits." -e 'my_hdr From:RFdiffusion (RFdiffusion)' -a top20_scoring_solutions.tar.bz2 -a animated.gif -- $USER@nih.gov
else
	echo -e "Files too large to attach to email."
	cat binders.html | mutt -e 'set content_type=text/html' -s "Binder design for $model_id is complete with $percent_hits% hits." -e 'my_hdr From:RFdiffusion (RFdiffusion)' -a animated.gif -- $USER@nih.gov
fi

# Clean up
rm chimovie*.png animated.gif
rm top20_scoring_solutions.tar.bz2
mkdir rfdiff_models
mkdir af2_models
mkdir silent_files
mkdir logs
mv rfdiff_*.pdb rfdiff_*.trb rfdiff_models/
mv rfdiff_models/rfdiff_*_af2pred.pdb af2_models/
mv *.silent* silent_files/
mv out.sc check.point *.txt *.cxc *.html *.log logs/

# Copy to Box, if account details are set up
if [ -e $HOME/.netrc ]; then
# Start lftp session
lftp ftp.box.com << EOF
cd RFdiffusion
mirror -R --no-symlinks $OUT_DIR/$project
bye
EOF
fi

echo -e "Done."

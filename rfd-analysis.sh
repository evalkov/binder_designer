#!/bin/bash
#SBATCH --job-name=bin_des_an
#SBATCH --nodes=1
#SBATCH --cpus-per-task=16
#SBATCH --mem-per-cpu=1G
#SBATCH --time=24:00:00
#SBATCH --mail-type=ALL
#SBATCH --mail-user=valkove2

model="/scratch/cluster_scratch/valkove2/rfdiffusion/caf1_20240223_181551/caf1.pdb"
PROC_DIR="/scratch/cluster_scratch/valkove2/rfdiffusion/caf1_20240223_181551"

cd $PROC_DIR

RFDIFFUSION_DIR="/home/valkove2/soft/RFdiffusion"

# Examine hydra config files for settings
hydra_config_1=`find . -type f -name config.yaml`
hydra_config_2=`find . -type f -name hydra.yaml`
hydra_config_3=`find . -type f -name overrides.yaml`
inference_log=`find . -type f -name run_inference.log` 


# Check if the variables are set (i.e., files were found)
if [ -z "$hydra_config_1" ] || [ -z "$hydra_config_2" ] || [ -z "$hydra_config_3" ] || [ -z "$inference_log" ]; then
    echo -e "Hydra config files are missing."
    exit 1
else
    echo -e "Found Hydra config files."
fi

device=`awk '/Found GPU with device_name/ {sub(/.*run RFdiffusion on /, ""); print}' $inference_log`
num_binders=`awk -F': ' '/num_designs/ {print $2}'  $hydra_config_1`
binder_length=`sed -n '/contigs:/,/inpaint_seq:/{/contigs:/d; /inpaint_seq:/d; p}' $hydra_config_1 | awk '{print $NF}'`
model_id=`basename $model`
hotspot_residue=`sed -n '/hotspot_res:/,/potentials/p' $hydra_config_1 | grep -o '[A-Za-z][0-9]\+' | paste -sd, -`
start_time=$(date -d "$(grep "output_dir:" $hydra_config_2 | awk -F/ '{print $NF}' | sed 's/-/:/g')" +"%I:%M %p")
start_date=$(date -d "$(grep "output_dir:" $hydra_config_2 | awk -F/ '{print $(NF-1)}')" +"%B %d, %Y")

if [ -e logs/out.sc ]; then
	finish_time=`date -d "$(stat -c '%y' "$(ls -lt logs/out.sc | head -n 1 | awk '{print $NF}')" | awk '{print $1 " " $2}')" +"%I:%M %p"`
	finish_date=`date -d "$(stat -c '%y' "$(ls -lt logs/out.sc | head -n 1 | awk '{print $NF}')" | awk '{print $1}')" +"%B %d, %Y"`
else
	finish_time=`date -d "$(stat -c '%y' "$(ls -lt out.sc | head -n 1 | awk '{print $NF}')" | awk '{print $1 " " $2}')" +"%I:%M %p"`
	finish_date=`date -d "$(stat -c '%y' "$(ls -lt out.sc | head -n 1 | awk '{print $NF}')" | awk '{print $1}')" +"%B %d, %Y"`
fi

echo "\
BINDER DESIGN PIPELINE
Written by Eugene Valkov, NCI/NIH.
" | tee -a $PROC_DIR/settings.txt

# Check if the hydra config file contains "ckpt_override_path"
if grep -q "ckpt_override_path" "$hydra_config_3"; then
    # If yes, print out the line containing that string
    betamodel=`grep "ckpt_override_path" $hydra_config_3 | awk -F/ '{print $NF}'`
    echo -e "\
RFdiffusion v `awk -F"[',]" '/version=/{print $2}' $RFDIFFUSION_DIR/setup.py` with Complex_beta_ckpt.pt model weight used to design" | tee -a $PROC_DIR/settings.txt
else
        echo -e "\
RFdiffusion v `awk -F"[',]" '/version=/{print $2}' $RFDIFFUSION_DIR/setup.py` with default model weight used to design" | tee -a $PROC_DIR/settings.txt
fi

echo "\
$num_binders binders of $binder_length residues targeting $hotspot_residue of $model_id
Sequence assigned with ProteinMPNN and validated with AlphaFold via dl_binder_design
Computed on `hostname` using $device GPU
All files at $PROC_DIR
Started on on $start_date at $start_time
Completed on $finish_date at $finish_time." | tee -a $PROC_DIR/settings.txt

# Check if the file out.sc exists in the current directory or in the logs directory
if [[ -e "out.sc" ]]; then
    input_file="out.sc"
    echo "File out.sc exists in the current directory."
elif [[ -e "logs/out.sc" ]]; then
    ln -s logs/out.sc out.sc
    input_file="out.sc"
    echo "File out.sc exists in the logs directory. Made a soft link"
else
    echo "File out.sc does not exist in both the current directory and logs."
    # Handle the case where the file doesn't exist, e.g., exit or set a default
fi


# Process the out.sc file
{
    # Extract and process the header: remove the first column
    head -n 1 "$input_file" | awk '{$1=""; print $0}' OFS="\t" | sed 's/^\t//'  | awk 'BEGIN {OFS="\t"} {print $1, $2, $3, $5, $10}'

    # Process the data: convert to tab-delimited, skip the first line, remove the first column,
    # filter for pae_interaction < 10, and sort by pae_interaction
    awk '{$1=$1; print}' OFS="\t" "$input_file" | cut -f2- | awk -F'\t' '$3 < 10' | sort -k3,3n | awk 'BEGIN {OFS="\t"} {print $1, $2, $3, $5, $10}'
} > pdb_list.txt


# Check if binders identified, and if not, then send an email and exit
if ! grep -q "_dldesign_0_cycle1_af2pred" pdb_list.txt; then
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

# Make a movie of top hits

echo "open $model_id" > movie.cxc

# Loop through each line of the input file, skipping the header line
awk 'NR > 1 && NR <= 11 {print $5 ".pdb"}' binders_list.txt | while read -r file; do
    # Check if the file exists
    if [[ -e $file ]]; then
            echo "open $file" >> movie.cxc
    elif [[ -e af2_models/$file ]]; then
	    echo "open af2_models/$file" >> movie.cxc
    fi
done

echo "\
set bgcolor white
cartoon style protein modeh tube rad 2 sides 24
cartoon style width 2 thick 0.2
rainbow chain palette RdYlBu-5
lighting flat shadows false intensity 0.3
matchmaker all to #1 pairing bs
del H
color #1 cornflower blue
show all cartoons
select clear
cofr #1
view cofr false
hide all models
show #1 models	
movie record size 1500,1500 format png supersample 4 directory $PROC_DIR
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
hide #$counter/B cartoons
show #$counter models
turn y 12 360 models #1-$counter
2dlab text '`echo $file | sed 's|.*/||; s|\.pdb$||'`' color black size 10 x .03 y .03
wait 30
stop
hide #$counter models
2dlab delete
" >> movie.cxc
    elif [[ -e af2_models/$file ]]; then
        echo "\
hide #$counter/B cartoons
show #$counter models
turn y 12 360 models #1-$counter
2dlab text '`echo af2_models/$file | sed 's|.*/||; s|\.pdb$||'`' color black size 10 x .03 y .03
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
hit_binders=$(tail -n +2 binders_list.txt | wc -l)

tar -cvjf top20_scoring_solutions.tar.bz2 top20_scoring_solutions/
tar -cvjf top_scoring_solutions.tar.bz2 top_scoring_solutions/

# Check files size before attachign to email
maxsize=$((7*1024*1024)) # 7 MB in bytes
size1=$(stat -c%s top20_scoring_solutions.tar.bz2)
size2=$(stat -c%s binders.html)
size3=$(stat -c%s animated.gif)
totalsize=$((size1 + size2 + size3))

# workaround with libcrypto not accessing the correct openssl
export LD_LIBRARY_PATH=/lib64:\$LD_LIBRARY_PATH

if [ $totalsize -lt $maxsize ]; then
	cat binders.html | mutt -e 'set content_type=text/html' -s "Binder design for $model_id is complete with $hit_binders hits." -e 'my_hdr From:RFdiffusion (RFdiffusion)' -a top20_scoring_solutions.tar.bz2 -a animated.gif -- $USER@nih.gov
else
	echo -e "Files too large to attach to email."
	cat binders.html | mutt -e 'set content_type=text/html' -s "Binder design for $model_id is complete with $hit_binders hits." -e 'my_hdr From:RFdiffusion (RFdiffusion)' -a animated.gif -- $USER@nih.gov
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
mirror -R --no-symlinks $PROC_DIR
bye
EOF
fi

echo -e "Done."

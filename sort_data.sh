#!/bin/bash

#SBATCH --nodes=1
#SBATCH --cpus-per-task=16
#SBATCH --mem-per-cpu=1G
#SBATCH --time=00:30:00
#SBATCH --mail-type=ALL
#SBATCH --mail-user=valkove2

cd not9_TESTING_DO_NOT_DELETE

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
cat filtered_sorted_output.txt

mkdir top_scoring_solutions

awk 'NR > 1 {print $NF ".pdb"}' filtered_sorted_output.txt | while IFS= read -r file; do
        cp "$file" top_scoring_solutions/
done

mkdir top20_scoring_solutions

# Initialize a counter
counter=0
awk 'NR > 1 {print $NF ".pdb"}' filtered_sorted_output.txt | while IFS= read -r file; do
# Increment the counter
    ((counter++))

    # Check if the counter has reached 20, if so, exit the loop
    if [[ $counter -eq 21 ]]; then
        break
    fi
    cp "$file" top20_scoring_solutions/
done

# Initialize a counter
counter=0
# Loop through each line of the input file, skipping the header line
while IFS=$'\t' read -r binder_aligned_rmsd pae_binder pae_interaction pae_target plddt_binder plddt_target plddt_total target_aligned_rmsd time description; do
    # Skip the header line
    if [[ $description == "description" ]]; then
        continue
    fi

    # Increment the counter
    ((counter++))

    # Check if the counter has reached 20, if so, exit the loop
    if [[ $counter -eq 21 ]]; then
        break
    fi

    # Check if the file exists
    filename="${description}.pdb"
    if [[ -e $filename ]]; then
        echo "open $filename" >> binders.cxc
    fi
done < filtered_sorted_output.txt

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


# Initialize a counter
counter=0

# Loop through each line of the input file, skipping the header line
while IFS=$'\t' read -r binder_aligned_rmsd pae_binder pae_interaction pae_target plddt_binder plddt_target plddt_total target_aligned_rmsd time description; do
    # Skip the header line
    if [[ $description == "description" ]]; then
        continue
    fi

    # Increment the counter
    ((counter++))

    # Check if the counter has reached 20, if so, exit the loop
    if [[ $counter -eq 21 ]]; then
        break
    fi

    # Check if the file exists
    filename="${description}.pdb"
    if [[ -e $filename ]]; then
	    echo "\
interfaces select #$counter/B contacting #$counter/A bothSides true
contacts #$counter/A restrict #$counter/B intraMol false
show sel atoms
select clear
" >> binders.cxc
    fi
done < filtered_sorted_output.txt

echo "\
delete H
color byhetero
" >> binders.cxc

cp binders.cxc top20_scoring_solutions/
mv binders.cxc top_scoring_solutions/


# Make a movie of top hits

counter=0
while IFS=$'\t' read -r binder_aligned_rmsd pae_binder pae_interaction pae_target plddt_binder plddt_target plddt_total target_aligned_rmsd time description; do
    if [[ $description == "description" ]]; then
        continue
    fi
    ((counter++))
    # Check if the counter has reached 6, if so, exit the loop
    if [[ $counter -eq 11 ]]; then
        break
    fi

    # Check if the file exists
    filename="${description}.pdb"
    if [[ -e $filename ]]; then
        echo "open $filename" >> movie.cxc
    fi
done < filtered_sorted_output.txt

echo "open not9.pdb" > movie.cxc

counter=0
while IFS=$'\t' read -r binder_aligned_rmsd pae_binder pae_interaction pae_target plddt_binder plddt_target plddt_total target_aligned_rmsd time description; do
    if [[ $description == "description" ]]; then
        continue
    fi
    ((counter++))
    # Check if the counter has reached 6, if so, exit the loop
    if [[ $counter -eq 11 ]]; then
        break
    fi

    # Check if the file exists
    filename="${description}.pdb"
    if [[ -e $filename ]]; then
        echo "\
open $filename" >> movie.cxc
    fi
done < filtered_sorted_output.txt

echo "\
set bgcolor white
cartoon style protein modeh tube rad 2 sides 24
cartoon style width 2 thick 0.2
rainbow chain palette RdYlBu-5
lighting flat shadows false intensity 0.3
matchmaker all to #1 pairing bs
select /B
cartoon hide sel
show #1 cartoons
color #1 cornflower blue
select #1:166
hide sel atoms
show sel cartoons
show sel atoms
style sel sphere
color sel byelement
select clear
cofr #1:166
view cofr false
hide all models
show #1 models
movie record size 1500,1500 format png supersample 4 directory /mnt/beegfs/valkov/rfdiffusion/not9_20231224_225546
" >> movie.cxc

counter=1
while IFS=$'\t' read -r binder_aligned_rmsd pae_binder pae_interaction pae_target plddt_binder plddt_target plddt_total target_aligned_rmsd time description; do
    if [[ $description == "description" ]]; then
        continue
    fi
    ((counter++))
    # Check if the counter has reached 6, if so, exit the loop
    if [[ $counter -eq 11 ]]; then
        break
    fi

    # Check if the file exists
    filename="${description}.pdb"
    if [[ -e $filename ]]; then
        echo "\
show #$counter models
turn y 12 360 models #1-$counter
2dlab text '$description' color black size 10 x .03 y .03
wait 30
stop
hide #$counter models
2dlab delete
" >> movie.cxc
    fi
done < filtered_sorted_output.txt

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

cp not9.pdb top20_scoring_solutions/
cp movie.cxc top20_scoring_solutions/

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

cp top_scoring_solutions.html top20_scoring_solutions/
cp top_scoring_solutions.html top_scoring_solutions/

tar -cvjf top20_scoring_solutions.tar.bz2 top20_scoring_solutions/
tar -cvjf top_scoring_solutions.tar.bz2 top_scoring_solutions/

# Check files size before attachign to email
maxsize=$((7*1024*1024)) # 7 MB in bytes
size=$(stat -c%s top20_scoring_solutions.tar.bz2)

size1=$(stat -c%s top20_scoring_solutions.tar.bz2)
size2=$(stat -c%s animated.gif)
totalsize=$((size1 + size2))

# workaround with libcrypto not accessing the correct openssl
export LD_LIBRARY_PATH=/lib64:\$LD_LIBRARY_PATH

if [ $totalsize -lt $maxsize ]; then
	echo -e "<img src=\"cid:animated.gif\" />" | mutt -e 'set content_type=text/html' -s "Binder design is complete with $percent_hits% hits." -e 'my_hdr From:RFdiffusion (RFdiffusion)' -a top20_scoring_solutions.tar.bz2 -a animated.gif -- $USER@nih.gov
else
	cat top_scoring_solutions.html | mutt -e 'set content_type=text/html' -s "Binder design is complete with $percent_hits% hits." -e 'my_hdr From:RFdiffusion (RFdiffusion)' -a animated.gif -- $USER@nih.gov
fi

# Clean up
rm chimovie*.png
#rm filtered_sorted_output.txt
#rm top20_scoring_solutions.tar.bz2

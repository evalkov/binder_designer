#!/bin/bash
#SBATCH --job-name=scoring
#SBATCH --cpus-per-task=48
#SBATCH --mem-per-cpu=8G
#SBATCH --time=00:30:00
#SBATCH --mail-type=ALL
#SBATCH --output=scoring_%j.out
#SBATCH --error=scoring_%j.err

# Exit on error, unset variables, and pipe failures
set -euo pipefail

# Load necessary modules
module load python

# Define the source directory
SOFT_DIR="/home/valkove2/soft/binder_designer"

# Define data directory and the AlphaFold score table
SOURCE_DATA_DIR="/mnt/beegfs/valkov/rfdiffusion/binder781-model3/binder781-model3-T135Wbba/af2pred_pdbs"
AF2_SCORE="/mnt/beegfs//valkov/rfdiffusion/binder781-model3/binder781-model3-T135Wbba/processing/af2_scores/out.sc"

# Define the Python interpreter path
PYTHON_BIN=$(which python3)

# Log the Python version being used (useful for debugging)
echo "Using Python interpreter at: $PYTHON_BIN"
$PYTHON_BIN --version

# Define a function to execute Python scripts with common prefixes
run_script() {
    $PYTHON_BIN "$SOFT_DIR/$1" "${@:2}"
}

# Execute the Python scripts sequentially using the function
run_script compute.py
run_script create_symlinks.py --source_dir "$SOURCE_DATA_DIR" --settings_file compute_settings.txt
run_script process_binders.py --af2_score "$AF2_SCORE" --settings_file compute_settings.txt
run_script pisa.py --settings_file compute_settings.txt
run_script merge_binders_contacts.py
run_script scaler.py
run_script heatmap_generation.py top_50_binders_weighted_minmax_scaler.csv
run_script heatmap_generation.py top_50_binders_weighted_standard_scaler.csv
run_script plots_generation.py binders_weighted_composite_scores_minmax_scaler.csv
run_script plots_generation.py binders_weighted_composite_scores_standard_scaler.csv
run_script top_binders.py top_50_binders_weighted_standard_scaler.csv top_50_binders_weighted_minmax_scaler.csv
run_script chimerax_script_generation.py

# Compress folders with predicitons and plots
tar -cvjf predictions.tar.bz2 predictions/
mkdir -p plots
mv *.eps plots/
tar -cvjf plots.tar.bz2 plots/

# Count lines excluding the first line in each file
hit_binders=$(tail -n +2 binders.csv | wc -l)

# Check files size before attaching to email
maxsize=$((7*1024*1024)) # 7 MB in bytes

# Calculate the combined size of both files
plotsize=$(stat -c%s plots.tar.bz2)
predsize=$(stat -c%s predictions.tar.bz2)
totalsize=$((plotsize + predsize))

# workaround with libcrypto not accessing the correct openssl
export LD_LIBRARY_PATH=/lib64:\$LD_LIBRARY_PATH

if [ $totalsize -lt $maxsize ]; then
        echo "Analysis complete." | mutt -e 'set content_type=text/html' -s "Identified $hit_binders binders." -e 'my_hdr From:RFdiffusion (RFdiffusion)' -a plots.tar.bz2 -a predictions.tar.bz2 -- $USER@nih.gov
else
        echo "Analysis complete." | mutt -e 'set content_type=text/html' -s "Identified $hit_binders binders." -e 'my_hdr From:RFdiffusion (RFdiffusion)' -- $USER@nih.gov
fi

# Tidy up
#rm *.pdb


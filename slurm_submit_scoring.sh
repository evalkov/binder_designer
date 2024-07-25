#!/bin/bash
#SBATCH --job-name=contacts
#SBATCH --nodes=1
#SBATCH --cpus-per-task=48
#SBATCH --mem-per-cpu=1G
#SBATCH --time=00:10:00
#SBATCH --mail-type=ALL
#SBATCH --mail-user=valkove2

module load python
module load ccp4/8.0.010

DIR=/home/valkove2/soft/binder_designer
SOFT=/mnt/

# This must be set for pisa to work and for pisa.cfg to be read
export CCP4=/mnt/nasapps/production/ccp4/8.0.010

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
check_file_exists "$DIR/top_50_binders_weighted_minmax_scaler.csv"
execute_command "`which python` $SOFT/heatmap_generation.py $DIR/top_50_binders_weighted_minmax_scaler.csv"

check_file_exists "$DIR/top_50_binders_weighted_standard_scaler.csv"
execute_command "`which python` $SOFT/heatmap_generation.py $DIR/top_50_binders_weighted_standard_scaler.csv"

check_file_exists "$DIR/binders_weighted_composite_scores_minmax_scaler.csv"
execute_command "`which python` $SOFT/plots_generation.py $DIR/binders_weighted_composite_scores_minmax_scaler.csv"

check_file_exists "$DIR/binders_weighted_composite_scores_standard_scaler.csv"
execute_command "`which python` $SOFT/plots_generation.py $DIR/binders_weighted_composite_scores_standard_scaler.csv"

# Check if both input files exist before running the compare_csv script
check_file_exists "$DIR/top_50_binders_weighted_standard_scaler.csv"
check_file_exists "$DIR/top_50_binders_weighted_minmax_scaler.csv"
execute_command "`which python` $SOFT/compare_csv.py $DIR/top_50_binders_weighted_standard_scaler.csv $DIR/top_50_binders_weighted_minmax_scaler.csv"

execute_command "`which python` $SOFT/top_binder_scatter_plot_generation.py"
execute_command "`which python` $SOFT/chimerax_script_generation.py"

# Tidy up
mkdir -p scoring_tables && mv *.csv scoring_tables/
mkdir plots && mv *.eps plots/

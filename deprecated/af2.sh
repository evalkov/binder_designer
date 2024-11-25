#!/bin/bash
#SBATCH --job-name=af2
#SBATCH --partition=gpu
#SBATCH --gres=gpu:v100:2
#SBATCH --nodes=1
#SBATCH --cpus-per-task=16
#SBATCH --mem-per-cpu=1G
#SBATCH --time=12:00:00
#SBATCH --mail-type=ALL

export DL_BINDER_DESIGN_DIR="/path/dl_binder_design"
export SILENT_TOOLS_DIR="/path/silent_tools"
export PROTEINMPNN_DIR="/path/ProteinMPNN"

module purge
export PATH=${SILENT_TOOLS_DIR}:$PATH
module load dl_binder_design/af2_binder_design
module load tensorRT/8.6.1-cuda12 cuda/12.0

# Check if RFDIR is set
if [ -z "$RFDIR" ]; then
  echo "RFDIR is not set. Exiting."
  exit 1
fi

cd $RFDIR

$DL_BINDER_DESIGN_DIR/af2_initial_guess/predict.py \
	-silent mpnn.silent \
	-outsilent af2.silent

silentextract af2.silent

# Export the RFDIR variable so it can be used in other scripts
export RFDIR

# Execute scoring.sh
/home/valkove2/scoring.sh


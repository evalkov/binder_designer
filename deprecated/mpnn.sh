#!/bin/bash
#SBATCH --job-name=mpnn
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

unset PYTHONPATH
module purge
export PATH=${SILENT_TOOLS_DIR}:$PATH
module load dl_binder_design/proteinmpnn_binder_design

rm check.point 2>/dev/null

# Check if RFDIR is set
if [ -z "$RFDIR" ]; then
  echo "RFDIR is not set. Exiting."
  exit 1
fi

cd $RFDIR

#add 'FIXED' labels to your pdbs which will be recognized by the ProteinMPNN scripts
$DL_BINDER_DESIGN_DIR/helper_scripts/addFIXEDlabels.py \
	--pdbdir pdb \
	--trbdir trb \
	--verbose

silentfrompdbs pdb/*.pdb > data.silent

$DL_BINDER_DESIGN_DIR/mpnn_fr/dl_interface_design.py \
	-silent data.silent \
	-relax_cycles 0 \
	-seqs_per_struct 2 \
	-outsilent mpnn.silent

sbatch --export=ALL,RFDIR="$RFDIR" --dependency=afterok:$SLURM_JOB_ID /home/valkove2/af2.sh

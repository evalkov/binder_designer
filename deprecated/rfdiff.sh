#!/bin/bash
#SBATCH --job-name=rfdiff
#SBATCH --partition=gpu
#SBATCH --gres=gpu:v100:2
#SBATCH --nodes=1
#SBATCH --cpus-per-task=16
#SBATCH --mem-per-cpu=1G
#SBATCH --time=12:00:00
#SBATCH --mail-type=ALL

export RFDIFFUSION_DIR="/path/RFdiffusion"

module purge
module load rfdiffusion/1.1.0

PROJECT_ID="binder781"
OUTPUT_DIR="/path"
TARGET_PDB="/path/model.pdb"

# Create a unique output directory
RANDOM_STRING=$(cat /dev/urandom | tr -dc 'a-zA-Z0-9' | fold -w 8 | head -n 1)
RFDIR="$OUTPUT_DIR/$PROJECT_ID-${RANDOM_STRING}"

mkdir -p $RFDIR
cd $RFDIR

# Define the parameters
echo "inference.output_prefix=$RFDIR/rfdiff_design \
inference.input_pdb=$TARGET_PDB \
'contigmap.contigs=[A1-44/9-15/0 C10-263]' \
'ppi.hotspot_res=[C40,C42,C230]' \
inference.num_designs=1000 \
inference.ckpt_override_path=$RFDIFFUSION_DIR/models/Complex_beta_ckpt.pt" > rfdiff.inp

# Run the script with the parameters
cat rfdiff.inp | xargs $RFDIFFUSION_DIR/scripts/run_inference.py

mkdir pdb && mv *.pdb pdb/
mkdir trb && mv *.trb trb/

sbatch --export=ALL,RFDIR="$RFDIR" --dependency=afterok:$SLURM_JOB_ID /home/valkove2/mpnn.sh

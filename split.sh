#!/bin/bash
#SBATCH --job-name=RFdiffusion
#SBATCH --partition=gpu
#SBATCH --gres=gpu:v100:2
#SBATCH --nodes=1
#SBATCH --cpus-per-task=16
#SBATCH --mem-per-cpu=1G
#SBATCH --time=7-00:00:00
#SBATCH --mail-type=ALL
#SBATCH --mail-user=valkove2

model="/home/valkove2/pdbs/not-module.pdb"
binder_length="50-100"
hotspot_residue="506"
num_binders="1000"
betamodel="no"

RFDIFFUSION_DIR="/home/valkove2/soft/RFdiffusion"
OUT_DIR="/scratch/cluster_scratch/valkove2/rfdiffusion"
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

silentfrompdbs rfdiff_*.pdb > data.silent

$DL_BINDER_DESIGN_DIR/mpnn_fr/dl_interface_design.py -silent data.silent -outsilent mpnn.silent | tee -a $OUT_DIR/$project/proteinmpnn.log

module purge
module load dl_binder_design/af2_binder_design
module load tensorRT/8.6.1-cuda12 cuda/12.0

$DL_BINDER_DESIGN_DIR/af2_initial_guess/predict.py -silent mpnn.silent -outsilent af2.silent | tee -a $OUT_DIR/$project/alphafold.log

silentextract af2.silent

cat out.sc | mutt -s "Binder design for $model_id complete." -e 'my_hdr From:RFdiffusion (RFdiffusion)' -- $USER@nih.gov

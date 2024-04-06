#!/bin/bash
#SBATCH --job-name=RFdiffusion
#SBATCH --partition=gpu
#SBATCH --gres=gpu:v100:2
#SBATCH --nodes=1
#SBATCH --cpus-per-task=16
#SBATCH --mem-per-cpu=1G
#SBATCH --time=4-00:00:00
#SBATCH --mail-type=ALL
#SBATCH --mail-user=valkove2

model="/home/valkove2/pabp.pdb"
binder_length="70-100"
hotspot_residue="116"
num_binders="1000"

RFDIFFUSION_DIR="/home/valkove2/soft/RFdiffusion"
OUT_DIR="/mnt/beegfs/valkov/rfdiffusion"
DL_BINDER_DESIGN_DIR="/home/valkove2/soft/dl_binder_design"
SILENT_TOOLS="/home/valkove2/soft/silent_tools"

model_id=`basename "$model" .pdb`
hotspot=`awk '/^ATOM/ && $6 == '$hotspot_residue' {residue_name = $4} END {print residue_name}' $model`
template_chain=`awk '/^ATOM/ {chain_id = substr($0, 22, 1)} END {print chain_id}' $model`
template_length=`awk '/^ATOM/ {residue = $6} END {print residue}' $model`
timestamp=$(date +"%Y%m%d_%H%M%S")
project="rfdiff_$timestamp"
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
	#inference.ckpt_override_path=$RFDIFFUSION_DIR/models/Complex_beta_ckpt.pt

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

echo -e "Generated binders are in $OUT_DIR/$project" | mutt -s "RFdiffusion for $model_id is finished" -e 'my_hdr From:RFdiffusion (RFdiffusion)' -b eugene.valkov@gmail.com -- "$USER"@nih.gov


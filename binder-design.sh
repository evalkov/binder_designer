#!/bin/bash
#SBATCH --job-name=bind_des
#SBATCH --partition=gpu
#SBATCH --gres=gpu:v100:2
#SBATCH --nodes=1
#SBATCH --cpus-per-task=16
#SBATCH --mem-per-cpu=1G
#SBATCH --time=3-00:00:00
#SBATCH --mail-type=ALL
#SBATCH --mail-user=$USER

PARAMS_FILE="$PWD/binder_designer.inp"

# edit these software locations

# this set of binder design and analysis scripts
BINDER_DESIGNER_DIR="$HOME/soft/binder_designer"

# below assumes installation in $HOME/soft/
RFDIFFUSION_DIR="$HOME/soft/RFdiffusion"
DL_BINDER_DESIGN_DIR="$HOME/soft/dl_binder_design"
SILENT_TOOLS="$HOME/soft/silent_tools"
PROTEINMPNN="$HOME/soft/ProteinMPNN"


# do not edit below

module purge
module load rfdiffusion/1.1.0

# Check if the input file with parameters exists
if [ -f "$PARAMS_FILE" ]; then
	# Extract values from parameters file
	OUT_DIR=$(grep 'output=' $PARAMS_FILE | cut -d'=' -f2)
	input_pdb=$(grep 'inference.input_pdb=' $PARAMS_FILE | cut -d'=' -f2)
	contigmap_contigs=$(grep 'contigmap.contigs=' $PARAMS_FILE | cut -d'=' -f2)
	ppi_hotspot_res=$(grep 'ppi.hotspot_res=' $PARAMS_FILE | cut -d'=' -f2)
    	num_designs=$(grep 'inference.num_designs=' $PARAMS_FILE | cut -d'=' -f2)

	timestamp=$(date +"%Y%m%d_%H%M%S")
	project=$(sed 's|.*/||; s|\.pdb$||' <<< "$input_pdb")_"$timestamp"

	mkdir $OUT_DIR/"$project"
	cp $input_pdb $OUT_DIR/$project
	cd $OUT_DIR/$project

	# Construct the command with values extracted from parameters file
    	cmd="$RFDIFFUSION_DIR/scripts/run_inference.py \
            'inference.output_prefix=$OUT_DIR/$project/rfdiff' \
            'inference.input_pdb=$input_pdb' \
            'contigmap.contigs=$contigmap_contigs' \
            'ppi.hotspot_res=$ppi_hotspot_res' \
            'inference.num_designs=$num_designs'"

	# Check if the use of beta model is specified in parameters file and conditionally append to cmd
	use_beta=$(grep 'Complex_beta_ckpt.pt' $PARAMS_FILE)
	if [ -n "$use_beta" ]; then
		cmd+=" 'inference.ckpt_override_path=$RFDIFFUSION_DIR/models/Complex_beta_ckpt.pt'"
	fi
	cp $PARAMS_FILE $OUT_DIR/$project
else
	echo "Error: Input file binder_designer.inp does not exist."
    	exit 1
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

RANDOM_STRING=$(cat /dev/urandom | tr -dc 'a-zA-Z0-9' | fold -w 8 | head -n 1)
ANALYSIS_SCR="binder-design-analysis.sh"
ANALYSIS_RUN="${ANALYSIS_SCR%.*}_$RANDOM_STRING.${ANALYSIS_SCR##*.}"

cp $BINDER_DESIGNER_DIR/binder-design-analysis.sh $ANALYSIS_RUN

awk -v out_dir="$OUT_DIR/$project" '/cd \$PROC_DIR/ {print "PROC_DIR=\""out_dir"\""; found=1} {print} END {if (!found) print "PROC_DIR=\""out_dir"\""}' $ANALYSIS_RUN > tmpfile && mv tmpfile $ANALYSIS_RUN
awk -v model="$input_pdb" '/cd \$PROC_DIR/ {print "model=\""model"\""; found=1} {print}' $ANALYSIS_RUN > tmpfile && mv tmpfile $ANALYSIS_RUN
awk -v rfdir="$RFDIFFUSION_DIR" '/cd \$PROC_DIR/ {print "RFDIFFUSION_DIR=\""rfdir"\""; found=1} {print}' $ANALYSIS_RUN > tmpfile && mv tmpfile $ANALYSIS_RUN

sbatch --dependency=afterok:$SLURM_JOB_ID $ANALYSIS_RUN


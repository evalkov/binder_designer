#!/bin/bash
#SBATCH --ntasks=1
#SBATCH --partition=gpu
#SBATCH --gres=gpu:v100:2
#SBATCH --nodes=1
#SBATCH --cpus-per-task=32
#SBATCH --open-mode=append
#SBATCH --time=3-00:00:00
#SBATCH --mem=200g
#SBATCH --mail-type=ALL
#SBATCH --mail-user=username

module load cudnn/8.8.3-cuda12
module load tensorRT/8.6.1-cuda12
module load evobind/v2

# Set the location of the target sequence file
fasta=/path/filename.fasta

# Create a unique output directory
RANDOM_STRING=$(cat /dev/urandom | tr -dc 'a-zA-Z0-9' | fold -w 8 | head -n 1)
OUTPUT_DIR="/path/evobind-${RANDOM_STRING}"
# Create output directory if it doesn't exist
mkdir -p $OUTPUT_DIR
# Redirect output and error to log.txt in OUTPUT_DIR
exec > >(tee -i ${OUTPUT_DIR}/log.txt)
exec 2>&1
# Copy the script itself into the OUTPUT_DIR
cp $0 $OUTPUT_DIR/
# Copy the target sequence into the OUTPUT_DIR
cp $fasta $OUTPUT_DIR/
cd $OUTPUT_DIR


#############PARAMETERS#############
BASE=/mnt/nasapps/production/evobind/v2/EvoBind
 #Change this depending on your local path
DATADIR=$OUTPUT_DIR #The output (designs) will also be written here
RECEPTORID=filename
###Receptor interface residues - provide with --receptor_if_residues=$RECEPTORIFRES if using
RECEPTORIFRES="40"
####Receptor fasta sequence####
RECEPTORFASTA=$DATADIR/$RECEPTORID'.fasta'
###Peptide length###
PEPTIDELENGTH=25
#NUMBER OF ITERATIONS
NITER=1000 #This will likely have to be modified depending on the outcome of the design

#########Step1: Create MSA with HHblits#########
HHBLITSDB=$BASE/data/uniclust30_2018_08/uniclust30_2018_08
#Write individual fasta files for all unique sequences
if test -f $DATADIR/$RECEPTORID'.a3m'; then
	echo $DATADIR/$RECEPTORID'.a3m' exists
else
	$BASE/hh-suite/build/bin/hhblits -i $RECEPTORFASTA -d $HHBLITSDB -E 0.001 -all -n 2 -oa3m $DATADIR/$RECEPTORID'.a3m'
fi
#MSA
MSA=$DATADIR/$RECEPTORID'.a3m'


#########Step2: Design binder#########
##### AF2 CONFIGURATION #####
PARAM=$BASE'/src/AF2/'
PRESET='full_dbs' #Choose preset model configuration - no ensembling (full_dbs) and (reduced_dbs) or 8 model ensemblings (casp14).
MAX_RECYCLES=8 #max_recycles (default=3)
MODEL_NAME='model_1' #model_1_ptm
MSAS="$MSA" #Comma separated list of msa paths

#Optimise a binder
#conda activate evobind
python3 $BASE/src/mc_design.py \
--receptor_fasta_path=$RECEPTORFASTA \
--peptide_length=$PEPTIDELENGTH \
--msas=$MSAS \
--output_dir=$DATADIR/ \
--model_names=$MODEL_NAME \
--data_dir=$PARAM \
--max_recycles=$MAX_RECYCLES \
--num_iterations=$NITER \
--predict_only=False \
--cyclic_offset=1
[valkove2@fsitgl-

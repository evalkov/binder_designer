#! /usr/bin/env bash
#SBATCH --partition=gpu
#SBATCH --gres=gpu:v100:2
#SBATCH --nodes=1
#SBATCH --cpus-per-task=16
#SBATCH --mem-per-cpu=20G
#SBATCH --time=12:00:00
#SBATCH --mail-type=ALL
#SBATCH --mail-user=valkove2

module purge
module load dl_binder_design/af2_binder_design
module load tensorRT/8.6.1-cuda12 cuda/12.0

mv out.silent mpnn.silent

/home/valkove2/soft/dl_binder_design/af2_initial_guess/predict.py -silent mpnn.silent

export PATH=/home/valkove2/soft/silent_tools:$PATH

silentextract out.silent

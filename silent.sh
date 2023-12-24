#! /usr/bin/env bash
#SBATCH --job-name=silent
#SBATCH --partition=gpu
#SBATCH --gres=gpu:v100:2
#SBATCH --nodes=1
#SBATCH --cpus-per-task=16
#SBATCH --mem-per-cpu=20G
#SBATCH --time=03:00:00
#SBATCH --mail-type=ALL
#SBATCH --mail-user=valkove2

export PATH=${HOME}/soft/silent_tools:$PATH

module purge
module load dl_binder_design/proteinmpnn_binder_design

silentfrompdbs /scratch/cluster_scratch/valkove2/rfdiff/rfdiff_20231223_230703/rfdiff_*.pdb > my_designs.silent

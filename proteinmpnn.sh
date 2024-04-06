#! /usr/bin/env bash
#SBATCH --job-name=protMPNN
#SBATCH --partition=gpu
#SBATCH --gres=gpu:v100:2
#SBATCH --nodes=1
#SBATCH --cpus-per-task=16
#SBATCH --mem-per-cpu=20G
#SBATCH --time=1-00:00:00
#SBATCH --mail-type=ALL
#SBATCH --mail-user=valkove2

export PATH=${HOME}/soft/silent_tools:$PATH

module purge
module load dl_binder_design/proteinmpnn_binder_design

/home/valkove2/soft/dl_binder_design/mpnn_fr/dl_interface_design.py -silent my_designs.silent

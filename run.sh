#!/bin/bash
#SBATCH --job-name=contacts
#SBATCH --nodes=1
#SBATCH --cpus-per-task=48
#SBATCH --mem-per-cpu=1G
#SBATCH --time=00:10:00
#SBATCH --mail-type=ALL
#SBATCH --mail-user=valkove2

module load python
module load ccp4/8.0.010

# This must be set for pisa to work and for pisa.cfg to be read
export CCP4=/mnt/nasapps/production/ccp4/8.0.010

`which python` /mnt/beegfs/valkov/rfdiffusion/binder781-model1/binder781-model1-TFldTo5K/parallel_pisa_estimates.py

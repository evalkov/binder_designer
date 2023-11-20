#!/bin/bash
#SBATCH --job-name=RFdiff
#SBATCH --output=rfdiff.out
#SBATCH --partition=gpu
#SBATCH --gres=gpu:v100:2
#SBATCH --nodes=1
#SBATCH --cpus-per-task=16
#SBATCH --mem-per-cpu=1G
#SBATCH --time=02:00:00
#SBATCH --mail-type=ALL
#SBATCH --mail-user=valkove2

module load rfdiffusion/1.1.0

# This will create a timestamp in the format of YYMMDD_HHMMSS
timestamp=$(date +"%Y%m%d_%H%M%S")

# Create a directory name with the timestamp
project="rfdiff_$timestamp"

model="/home/valkove2/caf1.pdb"
RFDIFFUSION_DIR="/home/valkove2/soft/RFdiffusion"
OUT_DIR="/scratch/cluster_scratch/valkove2/rfdiff"

# Create the directory
mkdir $OUT_DIR/"$project"

$RFDIFFUSION_DIR/scripts/run_inference.py \
	inference.output_prefix=$OUT_DIR/$project/ \
	inference.output_prefix=$OUT_DIR/$project/rfdiff \
	inference.input_pdb=$model \
	'contigmap.contigs=[C1-285/0 30-50]' \
       	'ppi.hotspot_res=[C40]' \
	inference.num_designs=10 \
	inference.ckpt_override_path=$RFDIFFUSION_DIR/models/Complex_beta_ckpt.pt

echo -e "RFDiffusion complete"

MOVIE_DIR="/mnt/projects/RNABL-GRID-SPE/active/Valkov/chimera_images"

echo "
set bgColor white
graphics rate maxFrameRate 30
movie record size 500,500 format png supersample 4 directory $MOVIE_DIR
open $OUT_DIR/$project/*.pdb
hide all models
~disp; cartoon style protein modeh tube rad 2 sides 24
rainbow chain palette RdYlBu-5
lighting flat shadows false intensity 1
" > $MOVIE_DIR/chimera_run.cxc

for file in $OUT_DIR/$project/*.pdb; do
        count=$((count + 1))
	echo "
show #$count models
wobble y 120 aspect 2 center cofr model #$count
2dlab text '#$count' color black size 25 x .03 y .03
wait 50
stop
hide #$count models
2dlab delete" >> $MOVIE_DIR/chimera_run.cxc
# Check if count is 5, then break the loop - reduces rendering time by ChimeraX
    if [ $count -eq 5 ]; then
        break
    fi
done

echo "
movie stop
show #1 models
save $MOVIE_DIR/session.cxs	
" >> $MOVIE_DIR/chimera_run.cxc

module load ChimeraX/1.5

ChimeraX --offscreen --script $MOVIE_DIR/chimera_run.cxc --exit

convert -dispose previous -delay 10 -loop 0 $MOVIE_DIR/chimovie*.png -coalesce $MOVIE_DIR/animated.gif

rm $MOVIE_DIR/chimovie*.png

# workaround with libcrypto not accessing the correct openssl 
export LD_LIBRARY_PATH=/lib64:$LD_LIBRARY_PATH

echo -e "<img src=\"cid:animated.gif\" />" | mutt -e 'set content_type=text/html' -s "Movie" -a $MOVIE_DIR/animated.gif -a $MOVIE_DIR/session.cxs -e 'my_hdr From:RFDiffusion (RFDiffusion)' -- "$USER"@nih.gov

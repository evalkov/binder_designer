#!/bin/bash
#SBATCH --job-name=RFdiffusion
#SBATCH --partition=gpu
#SBATCH --gres=gpu:v100:2
#SBATCH --nodes=1
#SBATCH --cpus-per-task=16
#SBATCH --mem-per-cpu=1G
#SBATCH --time=03:00:00
#SBATCH --mail-type=ALL
#SBATCH --mail-user=valkove2

module load rfdiffusion/1.1.0

model="/home/valkove2/not9.pdb"
binder_length="30-50"
hotspot_residue="166"
num_binders="10"


#### DO NOT EDIT BELOW ####

hotspot=`awk '/^ATOM/ && $6 == '$hotspot_residue' {residue_name = $4} END {print residue_name}' $model`
echo -e "Hostpot residue is $hotspot"

template_chain=`awk '/^ATOM/ {chain_id = substr($0, 22, 1)} END {print chain_id}' $model`

template_length=`awk '/^ATOM/ {residue = $6} END {print residue}' $model`

# This will create a timestamp in the format of YYMMDD_HHMMSS
timestamp=$(date +"%Y%m%d_%H%M%S")

# Create a directory name with the timestamp
project="rfdiff_$timestamp"

RFDIFFUSION_DIR="/home/valkove2/soft/RFdiffusion"
OUT_DIR="/scratch/cluster_scratch/valkove2/rfdiff"

# Create the directory
mkdir $OUT_DIR/"$project"

$RFDIFFUSION_DIR/scripts/run_inference.py \
	inference.output_prefix=$OUT_DIR/$project/rfdiff \
	inference.input_pdb=$model \
	'contigmap.contigs=['$template_chain''1'-'$template_length'/0 '$binder_length']' \
       	'ppi.hotspot_res=['$template_chain''$hotspot_residue']' \
	inference.num_designs=$num_binders \
	inference.ckpt_override_path=$RFDIFFUSION_DIR/models/Complex_beta_ckpt.pt

echo -e "RFDiffusion complete"

echo "
set bgcolor white
open $model
open $OUT_DIR/$project/rfdiff*.pdb
cartoon style protein modeh tube rad 2 sides 24
cartoon style width 2 thick 0.2
rainbow chain palette RdYlBu-5
lighting flat shadows false intensity 0.3
matchmaker all to #1 pairing bs
select /B
cartoon hide sel
show #1 cartoons
color #1 cornflower blue
sele #1:$hotspot_residue
cofr sel
hide sel atoms
show sel cartoons
show sel atoms
style sel sphere
color sel byhetero
select clear
hide all models
show #1 models
movie record size 1500,1500 format png supersample 4 directory $OUT_DIR/$project
" > $OUT_DIR/$project/chimera_run.cxc

count=1
for file in $OUT_DIR/$project/*.pdb; do
        count=$((count + 1))
	model_num=$((count - 1))
	echo "
show #$count models
turn y 12 360 models #1-$count
2dlab text '$model_num' color black size 25 x .03 y .03
wait 30
stop
hide #$count models
2dlab delete
" >> $OUT_DIR/$project/chimera_run.cxc
# Check if count is 11, then break the loop - reduces rendering time by ChimeraX
    if [ $count -eq 11 ]; then
        break
    fi
done

echo "
movie stop
save $OUT_DIR/$project/session.cxs	
" >> $OUT_DIR/$project/chimera_run.cxc

module load ChimeraX/1.5

ChimeraX --offscreen --script $OUT_DIR/$project/chimera_run.cxc --exit

convert \
	-dispose previous \
	-delay 10 \
       	-loop 0 \
       	-dither None \
	-colors 256 \
	-layers Optimize \
	-resize 350x350 \
	-filter Lanczos \
	-coalesce \
	$OUT_DIR/$project/chimovie*.png \
       	$OUT_DIR/$project/animated.gif

rm $OUT_DIR/$project/chimovie*.png

# workaround with libcrypto not accessing the correct openssl 
export LD_LIBRARY_PATH=/lib64:$LD_LIBRARY_PATH

echo -e "<img src=\"cid:animated.gif\" />" | mutt -e 'set content_type=text/html' -s "Movie" -a $OUT_DIR/$project/animated.gif -a $OUT_DIR/$project/session.cxs -e 'my_hdr From:RFDiffusion (RFDiffusion)' -b eugene.valkov@gmail.com -- "$USER"@nih.gov

if [ -e "$HOME"/.boxpassword ]; then
	username=`echo "$USER"@nih.gov`
	password=`awk '{print $1}' $HOME/.boxpassword`
	echo "\
set ftp:ssl-force true;
set mirror:parallel-directories true;
connect ftp://ftp.box.com;
user '$username' '$password';
cd RFdiffusion;
mirror -R --no-symlinks "$OUT_DIR"/"$project";
bye" > $HOME/."$project".lftpconfig
	chmod 600 $HOME/."$project".lftpconfig
	lftp -f $HOME/."$project".lftpconfig
	rm $HOME/."$project".lftpconfig
fi

echo -e "Transfer complete."

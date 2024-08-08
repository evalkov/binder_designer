#!/bin/bash

### GLOBAL PARAMS ###
export BINDER_DESIGNER_EXTERNAL="/mnt/beegfs/software/binder_designer/external/"
export CHUNKSIZE=50
export ARRAYSIZE=10
#####################

export RFDIFFUSION_DIR="$BINDER_DESIGNER_EXTERNAL/RFdiffusion"
export DL_BINDER_DESIGN_DIR="$BINDER_DESIGNER_EXTERNAL/dl_binder_design"
export SILENT_TOOLS_DIR="$BINDER_DESIGNER_EXTERNAL/silent_tools"
export PROTEINMPNN_DIR="$BINDER_DESIGNER_EXTERNAL/ProteinMPNN"

function print_centered {
    [[ $# == 0 ]] && return 1
    declare -i TERM_COLS="$(tput cols)"
    declare -i str_len="${#1}"
    [[ $str_len -ge $TERM_COLS ]] && {
        echo "$1";
        return 0;
    }
    declare -i filler_len="$(( (TERM_COLS - str_len) / 2 ))"
    [[ $# -ge 2 ]] && ch="${2:0:1}" || ch=" "
    filler=""
    for (( i = 0; i < filler_len; i++ )); do
        filler="${filler}${ch}"
    done
    printf "%s%s%s" "$filler" "$1" "$filler"
    [[ $(( (TERM_COLS - str_len) % 2 )) -ne 0 ]] && printf "%s" "${ch}"
    printf "\n"
    return 0
}

function print_header {
    echo
    print_centered "-" "-"
    print_centered "-" "-"
    print_centered "Binder Designer"
    print_centered "by E. Valkov"
    print_centered "and J. Caesar"
    print_centered "-" "-"
    print_centered "-" "-"
    echo
}

function print_footer {
    echo
    print_centered "-" "-"
    print_centered "-" "-"
    echo
}

function print_user_input {
    echo 
    echo "  Please see https://github.com/RosettaCommons/RFdiffusion for documentation about generating configmaps and selecting hotspot residues"
    echo
    echo -n "  Path to PDB file : "
    read pdbfile
    if test -f "$pdbfile"; then
        export PDBFILE=`readlink -f $pdbfile`
        export PROJECT_ID=`basename -- "$PDBFILE" ".pdb"`
        export OUTPUT_DIR="$(pwd)/$PROJECT_ID/$PROJECT_ID-$(cat /dev/urandom | tr -dc 'a-zA-Z0-9' | fold -w 8 | head -n 1)"
        # Get unique chain identifiers
        chains=$(grep '^ATOM' $PDBFILE | cut -c 22 | sort | uniq)
        echo
        # Iterate over each chain
        for chain in $chains; do
            # Get lines corresponding to the current chain
            chain_lines=$(grep '^ATOM' $PDBFILE | awk -v ch=$chain '$0 ~ " "ch" " {print}')
            # Get the first and last residue for the current chain
            first_residue=$(echo "$chain_lines" | head -n 1 | awk '{print $4, $6}')
            last_residue=$(echo "$chain_lines" | tail -n 1 | awk '{print $4, $6}')
            echo "    Chain $chain: First residue: $first_residue, Last residue: $last_residue"
        done
        echo
    else
        echo "ERROR : PDB file does not exist!"
        print_footer
        exit 1
    fi
    echo -n "  Contig map : "
    read contigmap
    if [ -z "${contigmap}" ]; then
        echo "ERROR : Contig map cannot be empty!"
        print_footer
        exit 1
    elif ! [[ "$contigmap" =~ ^[A-Z0-9\ \/-]+$ ]]; then
        echo "ERROR : Invalid characters in contig map!"
        print_footer
        exit 1    
    fi    
    export CONTIGMAP=${contigmap}
    export HOTSPOTRES
    while true; do
        echo -n "  Hotspot residue id [enter to exit] : "  
        read hotspotid
        if [ -z "${hotspotid}" ]; then
            break
        fi
        if [ -z "${HOTSPOTRES}" ]; then
            export HOTSPOTRES="${PDBTARGETCHAIN}$hotspotid"
        else
            export HOTSPOTRES="${HOTSPOTRES},${PDBTARGETCHAIN}$hotspotid"
        fi
    done
    while true; do
        echo -n "  RFDiffusion model (l to list) [Complex_beta] : "
        read rfdiffmodel
        if [ -z "${rfdiffmodel}" ]; then
            export RFDIFFMODEL="Complex_beta_ckpt.pt"
            break
        elif [ "${rfdiffmodel}" == "l" ]; then
            echo
            for file in $(ls $RFDIFFUSION_DIR/models/*_ckpt.pt); do
                echo "    $(basename -- $file '_ckpt.pt')"
            done
            echo
        elif test -f "$RFDIFFUSION_DIR/models/${rfdiffmodel}_ckpt.pt"; then
            export RFDIFFMODEL="${rfdiffmodel}_ckpt.pt"
            break
        else
            echo "ERROR : RFDiffusion model does not exist!"
            print_footer
            exit 1        
        fi
    done
    echo -n "  Number designs [1000]: "
    read ndesigns
    if [ -z "${ndesigns}" ]; then
        export NDESIGNS=1000
    elif [[ "${ndesigns}" =~ ^[0-9]+$ ]]; then
        export NDESIGNS=${ndesigns}
    else
        echo "ERROR: Number of designs is not a number"
        print_footer
        exit 1
    fi
    echo -n "  Number variants per design (Maximum 4) [1] : "
    read nvariants
    if [ -z "${nvariants}" ]; then
        export NVARIANTS=1
    elif [[ "${nvariants}" =~ ^[0-9]+$ ]]; then
        if [[ ${nvariants} -gt 4 ]]; then
            echo "ERROR: Number of variants should be 4 or less"
            print_footer
            exit 1       
        fi
        export NVARIANTS=${nvariants}
    else
        echo "ERROR: Number of variants is not a number"
        print_footer
        exit 1
    fi
}

function print_run_parameters {
    echo
    print_centered "-" "-"
    print_centered "run parameters"
    print_centered "-" "-"
    echo
    echo "  Project ID                 : $PROJECT_ID"
    echo "  Output directory           : $OUTPUT_DIR"
    echo "  PDB file                   : $PDBFILE"
    echo "  Contigmap                  : $CONTIGMAP"
    echo "  RFDiffusion model          : $RFDIFFMODEL"
    echo "  Hotspot residues           : $HOTSPOTRES"
    echo "  Number designs             : $NDESIGNS"
    echo "  Number variants per design : $NVARIANTS"
    echo
    echo "  Designs per chunk  : $CHUNKSIZE"
    echo "  Concurrent chunks  : $ARRAYSIZE"
    echo
    echo "  RFDiffusion directory       : $RFDIFFUSION_DIR"
    echo "  DL Binder Design directory  : $DL_BINDER_DESIGN_DIR"
    echo "  Silent Tools directory      : $SILENT_TOOLS_DIR"
    echo "  Protein MPNN directory      : $PROTEINMPNN_DIR"
}

function print_logline {
    dt=$(date '+%d/%m/%Y %H:%M:%S')
    echo "  $dt $1"
}

function print_logline_overwrite {
    dt=$(date '+%d/%m/%Y %H:%M:%S')
    echo -en "\r  $dt $1"
}

function write_sbatch_array_header {
    echo "#!/bin/bash"
    echo "#SBATCH --job-name=$1"
    echo "#SBATCH --partition=gpu"
    echo "#SBATCH --gres=gpu:v100:1"
    echo "#SBATCH --cpus-per-task=16"
    echo "#SBATCH --mem-per-cpu=1G"
    echo "#SBATCH --time=24:00:00"
    echo "#SBATCH --array=1-$2%$ARRAYSIZE"
    echo "#SBATCH --mail-type=FAIL,INVALID_DEPEND,TIME_LIMIT" 
}

function write_sbatch_header {
    echo "#!/bin/bash"
    echo "#SBATCH --job-name=$1"
    echo "#SBATCH --partition=gpu"
    echo "#SBATCH --gres=gpu:v100:1"
    echo "#SBATCH --cpus-per-task=16"
    echo "#SBATCH --mem-per-cpu=1G"
    echo "#SBATCH --time=24:00:00"
    echo "#SBATCH --mail-type=FAIL,INVALID_DEPEND,TIME_LIMIT" 
}

function generate_rfdiff_script {
    narray=$(($NDESIGNS / $CHUNKSIZE))
    remainder=$(($NDESIGNS % $CHUNKSIZE))
    ndesignarray=()
    if [[ $narray -gt 0 ]]; then
        for n in $(seq 1 $narray); do
            ndesignarray+=($CHUNKSIZE)
        done
    fi
    if [[ $remainder -gt 0 ]]; then
        ndesignarray+=($remainder)
    fi
    write_sbatch_array_header "rfdiff" ${#ndesignarray[*]}
    echo "NDESIGNARRAY=(${ndesignarray[*]})"
    echo "module purge"
    echo "module load rfdiffusion/1.1.0"
    echo 'echo "Starting RFDiffusion process ${SLURM_ARRAY_TASK_ID}. Please see rfdiff_${SLURM_ARRAY_TASK_ID}.log for detailed log" | tee -a binder_designer.log'
    echo "SECONDS=0"
    echo 'mkdir rfdiff_chunk_${SLURM_ARRAY_TASK_ID}'
    echo 'cd rfdiff_chunk_${SLURM_ARRAY_TASK_ID}'
    echo -n "$RFDIFFUSION_DIR/scripts/run_inference.py "
    echo -n "inference.output_prefix=$OUTPUT_DIR"'/rfdiff_chunk_${SLURM_ARRAY_TASK_ID}/rfdiff_design_${SLURM_ARRAY_TASK_ID} '
    echo -n "inference.input_pdb=$PDBFILE "
    echo -n "contigmap.contigs='[$CONTIGMAP]' "
    echo -n "ppi.hotspot_res='[$HOTSPOTRES]' "
    echo -n 'inference.num_designs=${NDESIGNARRAY[$(($SLURM_ARRAY_TASK_ID - 1))]} '
    echo -n "inference.ckpt_override_path=${RFDIFFUSION_DIR}/models/${RFDIFFMODEL} "
    echo '| tee -a '"${OUTPUT_DIR}"'/rfdiff_${SLURM_ARRAY_TASK_ID}.log'
    echo 'echo "RFDiffusion process ${SLURM_ARRAY_TASK_ID} finished" | tee -a '"$OUTPUT_DIR"'/binder_designer.log'
    echo 'echo "RFDiffusion process ${SLURM_ARRAY_TASK_ID} ran in ${SECONDS} seconds" | tee -a '"$OUTPUT_DIR"'/binder_designer.log'
    echo 'rm check.point 2>/dev/null'
    echo 'echo "Moving pdb and trb files for chunk ${SLURM_ARRAY_TASK_ID}" | tee -a '"$OUTPUT_DIR"'/binder_designer.log'
    echo 'mv *.pdb '"$OUTPUT_DIR"'/pdb/'
    echo 'mv *.trb $OUTPUT_DIR/trb/'
    echo 'if [[ ${SLURM_ARRAY_TASK_ID} -eq ${SLURM_ARRAY_TASK_MAX} ]]; then'
    echo "  cd $OUTPUT_DIR"
    echo '  sbatch --export=ALL --dependency=afterok:$SLURM_JOB_ID ./mpnn.sh'
    echo 'fi'
}

function generate_mpnn_script {
    write_sbatch_header "mpnn"
    echo 'unset PYTHONPATH'
    echo 'module purge'
    echo "export PATH=${SILENT_TOOLS_DIR}"':$PATH'
    echo 'module load dl_binder_design/proteinmpnn_binder_design'
    echo 'echo "Adding FIXED labels to pdbs for ProteinMPNN. Please see mpnn.log for detailed log" | tee -a binder_designer.log'
    echo "$DL_BINDER_DESIGN_DIR/helper_scripts/addFIXEDlabels.py --pdbdir pdb --trbdir trb --verbose | tee -a mpnn.log"
    echo 'echo "Converting pdbs to silent files for ProteinMPNN" | tee -a binder_designer.log'
    echo 'silentfrompdbs pdb/*.pdb > data.silent'
    echo 'echo "Starting ProteinMPNN process. Please see mpnn.log for detailed log" | tee -a binder_designer.log'
    echo 'SECONDS=0'
    echo "$DL_BINDER_DESIGN_DIR/mpnn_fr/dl_interface_design.py -silent data.silent -relax_cycles 0 -seqs_per_struct ${NVARIANTS} -outsilent mpnn.silent | tee -a mpnn.log"
    echo 'echo "ProteinMPNN finished" | tee -a binder_designer.log'
    echo 'echo "ProteinMPNN ran in ${SECONDS} seconds" | tee -a binder_designer.log'
    echo 'echo "Splitting mpnn.silent" | tee -a binder_designer.log'
    echo "silentsplit mpnn.silent $(($CHUNKSIZE * $NVARIANTS))"
    echo 'ls x* >> mpnn_chunks.txt'
    echo 'sbatch --export=ALL --dependency=afterok:$SLURM_JOB_ID ./af2.sh'
}

function generate_af2_script {
    narray=$(($NDESIGNS * $NVARIANTS / $CHUNKSIZE))
    remainder=$(($NDESIGNS * $NVARIANTS % $CHUNKSIZE))
    if [[ $remainder -gt 0 ]]; then
        write_sbatch_array_header "af2" $(($narray + 1))
    else
        write_sbatch_array_header "af2" $narray
    fi
    echo 'unset PYTHONPATH'
    echo 'module purge'
    echo "export PATH=${SILENT_TOOLS_DIR}"':$PATH'
    echo 'module load dl_binder_design/af2_binder_design'
    echo 'module load tensorRT/8.6.1-cuda12 cuda/12.0'
    echo 'echo "Starting af2 process ${SLURM_ARRAY_TASK_ID}. Please see af2_${SLURM_ARRAY_TASK_ID}.log for detailed log" | tee -a binder_designer.log'
    echo 'SECONDS=0'
    echo 'sed -n "${SLURM_ARRAY_TASK_ID}p" mpnn_chunks.txt | xargs '"$DL_BINDER_DESIGN_DIR"'/af2_initial_guess/predict.py  -checkpoint_name check_${SLURM_ARRAY_TASK_ID}.point -scorefilename out_${SLURM_ARRAY_TASK_ID}.sc -outsilent af2_${SLURM_ARRAY_TASK_ID}.silent -silent | tee -a '"${OUTPUT_DIR}"'/af2_${SLURM_ARRAY_TASK_ID}.log'
    echo 'echo "AF2 process ${SLURM_ARRAY_TASK_ID} finished" | tee -a binder_designer.log'
    echo 'echo "AF2 process ${SLURM_ARRAY_TASK_ID} ran in ${SECONDS} seconds" | tee -a binder_designer.log'
    echo 'echo "Extracting af2_${SLURM_ARRAY_TASK_ID}.silent" | tee -a binder_designer.log'
    echo 'silentextract af2_${SLURM_ARRAY_TASK_ID}.silent'
    echo 'if [[ ${SLURM_ARRAY_TASK_ID} -eq ${SLURM_ARRAY_TASK_MAX} ]]; then'
    echo "  cd $OUTPUT_DIR"
    echo '  sbatch --export=ALL --dependency=afterok:$SLURM_JOB_ID ./scoring.sh'
    echo 'fi'
}

print_header
print_user_input
print_run_parameters

echo
print_centered "-" "-"
print_centered "log"
print_centered "Press ctrl-c to exit. This will not kill running jobs"
print_centered "-" "-"
echo

print_logline "Creating directory structure"
mkdir -p $OUTPUT_DIR
mkdir -p $OUTPUT_DIR/pdb
mkdir -p $OUTPUT_DIR/trb

print_logline "Exporting script directory $(dirname -- $(readlink -f -- $0; ))"
export SCRIPT_DIR=$(dirname -- $(readlink -f -- $0; ))

print_logline "Moving to output directory"
cd $OUTPUT_DIR
touch binder_designer.log

print_logline "Saving run parameters"
print_run_parameters > run_parameters.txt
echo "<html><body><pre>" > run_parameters.html
print_run_parameters >> run_parameters.html
echo "</pre></body></html>" >> run_parameters.html

print_logline "Generating RFDiffusion script"
generate_rfdiff_script > rfdiff.sh

print_logline "Generating ProteinMPNN script"
generate_mpnn_script > mpnn.sh

print_logline "Generating AF2 script"
generate_af2_script > af2.sh

print_logline "Copying scoring script"
cp $SCRIPT_DIR/scoring.sh scoring.sh

print_logline "Submitting first job in the pipeline"
sbatch --export=ALL ./rfdiff.sh

CURRENTPROCESS="submission"
RFDIFFNPROC=0
AF2NPROC=0

while true; do

    if [ "$CURRENTPROCESS" == "rfdiff" ]; then
        nprocessed=$(ls $OUTPUT_DIR/rfdiff_chunk_* | grep pdb | wc -l)
        if(( $nprocessed > $RFDIFFNPROC )); then
            RFDIFFNPROC=$nprocessed
            print_logline_overwrite "RFDiffusion generated  $RFDIFFNPROC / $NDESIGNS designs"
        fi
   # elif [ "$CURRENTPROCESS" != "mpnn" ]; then
    #    print_logline_overwrite "MPNN generated $(ls $OUTPUT_DIR/rfdiff_chunk_*/*.pdb | wc -l) / $(( $NDESIGNS * $NVARIANTS )) design variants"
    elif [ "$CURRENTPROCESS" == "af2" ]; then
        nprocessed=0
        for file in $(ls $OUTPUT_DIR | grep out_ | grep .sc); do
            nprocessed=$(( $nprocessed + $(tail -n +2 $file | wc -l) ))
        done
        if(( $nprocessed > $AF2NPROC )); then
            AF2NPROC=$nprocessed
            print_logline_overwrite "Alphafold 2 processed ${AF2NPROC} / $(( $NDESIGNS * $NVARIANTS )) design variants"
        fi
    fi
    
    if test -f "$OUTPUT_DIR/predictions.tar.bz2"; then
        echo
        print_logline "Scoring process complete"
        break
    elif test -f "$OUTPUT_DIR/out.sc"; then
        if [ "$CURRENTPROCESS" != "scoring" ]; then
            CURRENTPROCESS=scoring
            echo
            print_logline "Scoring process started"
        fi
    elif test -f "$OUTPUT_DIR/out_1.sc"; then
        if [ "$CURRENTPROCESS" != "af2" ]; then
            CURRENTPROCESS=af2
            echo
            print_logline "Alphafold 2 process started. Please see af2_*.log for detailed logs"
        fi
    elif test -f "$OUTPUT_DIR/data.silent"; then
        if [ "$CURRENTPROCESS" != "mpnn" ]; then
            CURRENTPROCESS=mpnn
            echo
            print_logline "MPNN process started. Please see mpnn.log for detailed log"
        fi
    elif test -d "$OUTPUT_DIR/rfdiff_chunk_1"; then
        if [ "$CURRENTPROCESS" != "rfdiff" ]; then
            CURRENTPROCESS=rfdiff
            echo
            print_logline "RFDiffusion process started. Please see rfdiff_*.log for detailed logs"
        fi
    fi
    sleep 60
done

print_footer

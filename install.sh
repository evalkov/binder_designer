#!/bin/bash

# Set the installation directory
export BINDER_DESIGNER_DIR="${PWD}"
INSTALL_DIR="$PWD/soft"

# Ensure the installation directory exists
mkdir -p $INSTALL_DIR

# Navigate to the installation directory
cd $INSTALL_DIR

# Clone repositories
git clone https://github.com/RosettaCommons/RFdiffusion.git
git clone https://github.com/nrbennet/dl_binder_design.git
git clone https://github.com/bcov77/silent_tools.git
git clone https://github.com/dauparas/ProteinMPNN.git

# Content to be inserted
read -r -d '' PATHS << EOF
export BINDER_DESIGNER_DIR="${BINDER_DESIGNER_DIR}"
export RFDIFFUSION_DIR="${INSTALL_DIR}/RFdiffusion"
export DL_BINDER_DESIGN_DIR="${INSTALL_DIR}/dl_binder_design"
export SILENT_TOOLS_DIR="${INSTALL_DIR}/silent_tools"
export PROTEINMPNN_DIR="${INSTALL_DIR}/ProteinMPNN"
EOF

# Use awk to insert the content before "module purge"
awk -v insertion="$PATHS" '
    /module purge/ && !inserted {
        print insertion;
        inserted=1;
    }
    { print }
' "${BINDER_DESIGNER_DIR}/binder-design.sh" > temp && mv temp "${BINDER_DESIGNER_DIR}/binder-design.sh"

chmod +x "${BINDER_DESIGNER_DIR}/binder-design.sh"

# Use awk to insert the content before "cd $PROC_DIR"
awk -v insertion="export RFDIFFUSION_DIR=\"${INSTALL_DIR}/RFdiffusion\"" '
    /cd \$PROC_DIR/ && !inserted {
        print insertion;
        inserted=1;
    }
    { print }
' "${BINDER_DESIGNER_DIR}/binder-design-analysis.sh" > temp && mv temp "${BINDER_DESIGNER_DIR}/binder-design-analysis.sh"

chmod +x "${BINDER_DESIGNER_DIR}/binder-design-analysis.sh"

# Replace the placeholder with the full path of RFdiffusion installed above
sed -i "s|schedule_directory_path: null|schedule_directory_path: ${INSTALL_DIR}/RFdiffusion|" ${INSTALL_DIR}/RFdiffusion/config/inference/base.yaml

# Make a soft link to ProteinMPNN within the dl_binder_design directory
ln -s ${INSTALL_DIR}/ProteinMPNN ${INSTALL_DIR}/dl_binder_design/mpnn_fr/ProteinMPNN

# Download and extract AlphaFold2 model weights
mkdir -p ${INSTALL_DIR}/dl_binder_design/af2_initial_guess/model_weights/params
cd ${INSTALL_DIR}/dl_binder_design/af2_initial_guess/model_weights/params
wget https://storage.googleapis.com/alphafold/alphafold_params_2022-12-06.tar
tar --extract --verbose --file=alphafold_params_2022-12-06.tar

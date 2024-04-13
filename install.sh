#!/bin/bash

# Set the installation directory
INSTALL_DIR="$HOME/testsoft"

# Ensure the installation directory exists
mkdir -p $INSTALL_DIR

export BINDER_DESIGNER_DIR="${PWD}"

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

# Download the model weights into the RFDiffusion directory
mkdir -p ${INSTALL_DIR}/RFdiffusion/models
cd ${INSTALL_DIR}/RFdiffusion/models
wget http://files.ipd.uw.edu/pub/RFdiffusion/6f5902ac237024bdd0c176cb93063dc4/Base_ckpt.pt
wget http://files.ipd.uw.edu/pub/RFdiffusion/e29311f6f1bf1af907f9ef9f44b8328b/Complex_base_ckpt.pt
wget http://files.ipd.uw.edu/pub/RFdiffusion/60f09a193fb5e5ccdc4980417708dbab/Complex_Fold_base_ckpt.pt
wget http://files.ipd.uw.edu/pub/RFdiffusion/74f51cfb8b440f50d70878e05361d8f0/InpaintSeq_ckpt.pt
wget http://files.ipd.uw.edu/pub/RFdiffusion/76d00716416567174cdb7ca96e208296/InpaintSeq_Fold_ckpt.pt
wget http://files.ipd.uw.edu/pub/RFdiffusion/5532d2e1f3a4738decd58b19d633b3c3/ActiveSite_ckpt.pt
wget http://files.ipd.uw.edu/pub/RFdiffusion/12fc204edeae5b57713c5ad7dcb97d39/Base_epoch8_ckpt.pt
wget http://files.ipd.uw.edu/pub/RFdiffusion/f572d396fae9206628714fb2ce00f72e/Complex_beta_ckpt.pt
wget http://files.ipd.uw.edu/pub/RFdiffusion/1befcb9b28e2f778f53d47f18b7597fa/RF_structure_prediction_weights.pt

# Download and extract AlphaFold2 model weights
mkdir -p ${INSTALL_DIR}/dl_binder_design/af2_initial_guess/model_weights/params
cd ${INSTALL_DIR}/dl_binder_design/af2_initial_guess/model_weights/params
wget https://storage.googleapis.com/alphafold/alphafold_params_2022-12-06.tar
tar --extract --verbose --file=alphafold_params_2022-12-06.tar

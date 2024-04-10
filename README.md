The script expects an input file in the path, binder_designer.inp, which has the following content

output=/path-to-output
inference.input_pdb=/full-path-to-pdb-file
contigmap.contigs=[B1-268/0 20-30]
ppi.hotspot_res=[B166]
inference.num_designs=20

An additional line may specify a "beta" model, which generates a greater diversity of topologies.
Complex_beta_ckpt.pt

contigmap.contigs=[B1-268/0 20-30] specifies that the target model has 1-268 residues in chain B, "/0" specifies chain breaks, and 20-30 is the residue range for binders.

ppi.hotspot_res=[B166] where where B is the chain ID in the input pdb file of the hotspot residue and 166 is the residue number in the input pdb file of the hotspot residue.

Multiple hotspot residues can be specified as [A30,A33,A34], for example.

inference.num_designs=20 specifies 20 binder designs to generate.



INSTALLATION

Below assumes installation in $HOME/soft/

cd ~/soft
git clone https://github.com/RosettaCommons/RFdiffusion.git
git clone https://github.com/nrbennet/dl_binder_design.git
git clone https://github.com/bcov77/silent_tools.git
git clone https://github.com/dauparas/ProteinMPNN.git

Replace the /home/valkove2/soft/RFdiffusion with the full path of RFdiffusion installed above.

sed -i "s|schedule_directory_path: null|schedule_directory_path: \/home/valkove2/soft/RFdiffusion|" ~/soft/RFdiffusion/config/inference/base.yaml

Make a soft link to ProteinMPNN within the dl_binder_design directory.

ln -s ~/soft/ProteinMPNN ~/soft/dl_binder_design/mpnn_fr/ProteinMPNN

Download and extract AlphaFold2 model weights.

mkdir -p ~/soft/dl_binder_design/af2_initial_guess/model_weights/params && cd ~/soft/dl_binder_design/af2_initial_guess/model_weights/params
wget https://storage.googleapis.com/alphafold/alphafold_params_2022-12-06.tar
tar --extract --verbose --file=alphafold_params_2022-12-06.tar




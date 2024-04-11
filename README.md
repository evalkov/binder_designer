The script expects an input file in the path, binder_designer.inp, which has the following content

output=/path-to-output

inference.input_pdb=/full-path-to-pdb-file

contigmap.contigs=[B1-268/0 20-30]

ppi.hotspot_res=[B166]

inference.num_designs=20

An additional line may specify a "beta" model, which generates a greater diversity of topologies.

Complex_beta_ckpt.pt

contigmap.contigs=[B1-268/0 20-30] specifies that the target model has 1-268 residues in chain B, "/0" specifies chain breaks, and 20-30 is the residue range for binders.

ppi.hotspot_res=[B166], where B is the chain ID in the hotspot residue's input pdb file and 166 is the residue number.

Multiple hotspot residues can be specified as [A30,A33,A34], for example.

inference.num_designs=20 specifies 20 binder designs to generate.



INSTALLATION

Edit install.sh and specify the install directory for the repositories above. By default, they will be installed in ~/soft

then 

chmod +x install.sh

./install.sh


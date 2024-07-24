import os
import pandas as pd
import tarfile

# Load the dataset
file_path = 'top50_common.csv'
data = pd.read_csv(file_path)

# Specify the paths
predictions_dir = 'predictions'
binders_cxc_path = os.path.join(predictions_dir, 'binders.cxc')

# Create the predictions directory if it doesn't exist
os.makedirs(predictions_dir, exist_ok=True)

# Create the binders.cxc file inside the predictions directory and write the initial set of commands
with open(binders_cxc_path, 'w') as file:
    file.write("set bgcolor white\n")

# Process all binders
for idx, row in data.iterrows():
    binder_file = row['description'] + ".pdb"
    
    # Check if the file exists in the current directory or af2_models directory
    if os.path.exists(binder_file):
        with open(binders_cxc_path, 'a') as file:
            file.write(f"open {binder_file}\n")
        os.system(f"cp {binder_file} {predictions_dir}/")
    elif os.path.exists(f"af2_models/{binder_file}"):
        with open(binders_cxc_path, 'a') as file:
            file.write(f"open af2_models/{binder_file}\n")
        os.system(f"cp af2_models/{binder_file} {predictions_dir}/")

# Add the remaining commands to the binders.cxc file
with open(binders_cxc_path, 'a') as file:
    file.write("""cartoon style protein modeh tube rad 2 sides 24
cartoon style width 2 thick 0.2
rainbow chain palette RdYlBu-5
lighting simple shadows false intensity 0.5
view all
hide atoms
show cartoons
hide all models
show #1 models
matchmaker all to #1/B pairing bs
""")

# Initialize a counter
counter = 0

# Loop through all binders again to add interface and contact commands
for idx, row in data.iterrows():
    binder_file = row['description'] + ".pdb"
    counter += 1

    # Check if the file exists in the current directory or af2_models directory
    if os.path.exists(binder_file) or os.path.exists(f"af2_models/{binder_file}"):
        with open(binders_cxc_path, 'a') as file:
            file.write(f"""interfaces select #{counter}/B contacting #{counter}/A bothSides true
contacts #{counter}/A restrict #{counter}/B intraMol false
show sel atoms
select clear
""")

# Add the final commands to the binders.cxc file
with open(binders_cxc_path, 'a') as file:
    file.write("""delete H
color byhetero
""")

# Compress the predictions folder as tar.bz2
with tarfile.open('predictions.tar.bz2', 'w:bz2') as tar:
    tar.add(predictions_dir, arcname=os.path.basename(predictions_dir))

print(f"binders.cxc file has been created and populated successfully at {binders_cxc_path}.")
print(f"The directory {predictions_dir} has been compressed into predictions.tar.bz2.")

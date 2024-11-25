import pandas as pd

# Input files
salt_bridges_file = "qualifying_salt_bridges.csv"
out_sc_file = "out.sc"
output_file = "filtered_out.sc"

# Step 1: Load PDB_IDs from qualifying_salt_bridges.csv
print("Reading qualifying_salt_bridges.csv...")
salt_bridges_df = pd.read_csv(salt_bridges_file)

# Extract unique PDB_IDs and ensure they are clean
pdb_ids = set(salt_bridges_df["PDB_ID"].str.strip())
print(f"Extracted {len(pdb_ids)} unique PDB_IDs.")

# Step 2: Read and process the out.sc file
print("Reading out.sc...")
with open(out_sc_file, "r") as f:
    lines = f.readlines()

# Separate header and data lines
header_line = lines[0]  # The first line is the header
data_lines = lines[1:]  # The remaining lines are data

print(f"Found 1 header line and {len(data_lines)} data lines.")

# Step 3: Filter data_lines based on matching description field
filtered_lines = []
for line in data_lines:
    # Split line into fields and get the last field as description
    parts = line.strip().split()
    if len(parts) > 0:
        description = parts[-1].strip()  # Last field is the description
        if description in pdb_ids:
            filtered_lines.append(line)  # Add line if description matches PDB_ID

print(f"Filtered {len(filtered_lines)} matching lines.")

# Step 4: Write the filtered output
print("Writing filtered output...")
with open(output_file, "w") as f:
    f.write(header_line)  # Write the header
    f.writelines(filtered_lines)  # Write filtered data

print(f"Filtered out.sc file saved to {output_file}")


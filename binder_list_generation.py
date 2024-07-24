import pandas as pd
import os

# Get the current working directory
cwd = os.getcwd()

# Construct the file paths
contacts_path = os.path.join(cwd, 'contacts.csv')
binders_list_path = os.path.join(cwd, 'binders_list.txt')
final_binders_list_path = os.path.join(cwd, 'final_binders_list.csv')

# Read the CSV file
contacts_df = pd.read_csv(contacts_path)

# Read the TSV file
binders_list_df = pd.read_csv(binders_list_path, sep='\t')

# Combine the dataframes on the 'binder' column from contacts_df and 'description' column from binders_list_df
combined_df = pd.merge(binders_list_df, contacts_df[['binder', 'bsa_score', 'salt_bridges', 'h_bonds']],
                       left_on='description', right_on='binder', how='inner')

# Drop the duplicate 'binder' column after merging
combined_df.drop(columns=['binder'], inplace=True)

# Save the combined dataframe to a new CSV file
combined_df.to_csv(final_binders_list_path, index=False)

print("The combined CSV file has been created successfully in the current working directory.")


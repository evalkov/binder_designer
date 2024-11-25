import pandas as pd

# Load the CSV files
df1 = pd.read_csv('top_50_binders_weighted_standard_scaler.csv')
df2 = pd.read_csv('top_50_binders_weighted_minmax_scaler.csv')

# Merge the dataframes on the 'binder_seq' column
merged_df = pd.merge(df1, df2, on='binder_seq', suffixes=('_standard', '_minmax'))

# Select relevant columns and rename them to match the example output
merged_df = merged_df[['description_standard', 'binder_seq', 'weighted_composite_score_standard', 'weighted_composite_score_minmax']]
merged_df.columns = ['description', 'binder_seq', 'standard_scale_score', 'minmax_scale_score']

# Save the merged dataframe to a new CSV file
merged_df.to_csv('top50_common.csv', index=False)

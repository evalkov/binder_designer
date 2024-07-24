import pandas as pd
from sklearn.preprocessing import StandardScaler, MinMaxScaler

# Function to filter data based on residue and range input
def filter_data_by_residue_and_range(data, residue, start_range, end_range):
    filtered_data = data[data['binder_seq'].str[start_range-1:end_range].str.contains(residue)]
    return filtered_data

# Load the dataset
file_path = 'final_binders_list.csv'
data = pd.read_csv(file_path)

# Set the residue and range
residue_input = 'P/44-48'  # Example residue and range

# Parse the residue and range
residue, range_input = residue_input.split('/')
start_range, end_range = map(int, range_input.split('-'))
if start_range > end_range:
    raise ValueError("Start range cannot be greater than end range.")
data = filter_data_by_residue_and_range(data, residue, start_range, end_range)

if data.empty:
    raise ValueError("No entries match the specified residue and range.")

# Define weights for each metric
weights = {
    'binder_aligned_rmsd': 1,
    'pae_binder': 1,
    'pae_interaction': 1,
    'plddt_binder': 1,
    'bsa_score': 3,
    'salt_bridges': 3,
    'h_bonds': 1
}

# Identify metrics to reverse (lower is better)
reverse_metrics = ['binder_aligned_rmsd', 'pae_binder', 'pae_interaction']
# Identify metrics to normalize (higher is better)
higher_better_metrics = ['plddt_binder', 'bsa_score', 'salt_bridges', 'h_bonds']

# Reverse the metrics where lower is better
data[reverse_metrics] = data[reverse_metrics].apply(lambda x: -x)

# Combine the metrics for normalization
metrics_to_normalize = reverse_metrics + higher_better_metrics

# Apply weights
for metric in metrics_to_normalize:
    data[metric] *= weights[metric]

# Standard Scaler
standard_scaler = StandardScaler()
data_standard_scaled = data.copy()
data_standard_scaled[metrics_to_normalize] = standard_scaler.fit_transform(data[metrics_to_normalize])
data_standard_scaled[metrics_to_normalize] = data_standard_scaled[metrics_to_normalize].round(2)
data_standard_scaled['weighted_composite_score'] = data_standard_scaled[metrics_to_normalize].mean(axis=1)
data_standard_scaled.to_csv('binders_weighted_composite_scores_standard_scaler.csv', index=False, float_format='%.2f')

# MinMax Scaler
minmax_scaler = MinMaxScaler()
data_minmax_scaled = data.copy()
data_minmax_scaled[metrics_to_normalize] = minmax_scaler.fit_transform(data[metrics_to_normalize])
data_minmax_scaled[metrics_to_normalize] = data_minmax_scaled[metrics_to_normalize].round(2)
data_minmax_scaled['weighted_composite_score'] = data_minmax_scaled[metrics_to_normalize].mean(axis=1)
data_minmax_scaled.to_csv('binders_weighted_composite_scores_minmax_scaler.csv', index=False, float_format='%.2f')

# Select top 50 based on the new weighted composite score
top_50_standard_scaled = data_standard_scaled.nlargest(50, 'weighted_composite_score')
top_50_standard_scaled.to_csv('top_50_binders_weighted_standard_scaler.csv', index=False, float_format='%.2f')

top_50_minmax_scaled = data_minmax_scaled.nlargest(50, 'weighted_composite_score')
top_50_minmax_scaled.to_csv('top_50_binders_weighted_minmax_scaler.csv', index=False, float_format='%.2f')

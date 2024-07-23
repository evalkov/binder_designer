import pandas as pd
from sklearn.preprocessing import MinMaxScaler
import matplotlib.pyplot as plt
import seaborn as sns
import matplotlib.font_manager as fm

# Function to filter data based on residue and range input
def filter_data_by_residue_and_range(data, residue, start_range, end_range):
    filtered_data = data[data['binder_seq'].str[start_range-1:end_range].str.contains(residue)]
    return filtered_data

# Function to process y-axis labels to show only the differing text
def get_unique_suffixes(labels):
    # Split each label into parts
    split_labels = [label.split('_') for label in labels]
    # Transpose the list to get columns
    transposed = list(zip(*split_labels))
    
    unique_suffixes = []
    for parts in transposed:
        if len(set(parts)) > 1:
            unique_suffixes.append(parts)
        else:
            unique_suffixes.append([''] * len(parts))
    
    # Transpose back and join parts to form new labels
    processed_labels = ['_'.join(filter(None, parts)) for parts in zip(*unique_suffixes)]
    return processed_labels

# Load the dataset
file_path = 'final_binders_list.csv'
data = pd.read_csv(file_path)

# Ask user for specific residue and range
while True:
    residue_input = input("Enter the residue and range (e.g., W/20-25) or press Enter to skip: ")
    if residue_input:
        try:
            residue, range_input = residue_input.split('/')
            start_range, end_range = map(int, range_input.split('-'))
            if start_range > end_range:
                raise ValueError("Start range cannot be greater than end range.")
            data = filter_data_by_residue_and_range(data, residue, start_range, end_range)
            
            if data.empty:
                raise ValueError("No entries match the specified residue and range.")
            break
        except ValueError as e:
            print(f"Error: {e}. Please specify a residue and a valid range.")
    else:
        break

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

# Apply weights and then normalize
for metric in metrics_to_normalize:
    data[metric] *= weights[metric]
    
# Normalize the weighted metrics to keep them in the 0-1 scale
scaler = MinMaxScaler()
data[metrics_to_normalize] = scaler.fit_transform(data[metrics_to_normalize])

# Recalculate the weighted composite score
data['weighted_composite_score'] = data[metrics_to_normalize].mean(axis=1)

# Save the complete dataframe with composite scores to a CSV file, rounded to 3 decimal places
output_csv_path = 'top50_weighted_composite_scores.csv'
data.to_csv(output_csv_path, index=False, float_format='%.3f')

# Select top 50 based on the new weighted composite score
top_50_weighted_data = data.nlargest(50, 'weighted_composite_score')

# Set the description as the index for better readability
top_50_weighted_data.set_index('description', inplace=True)

# Process y-axis labels to show only the differing text
processed_labels = get_unique_suffixes(top_50_weighted_data.index)

# Define a custom color map similar to the provided image
custom_cmap = sns.color_palette("Blues", as_cmap=True)

# Select the columns used in normalization for the heatmap
heatmap_weighted_data = top_50_weighted_data[metrics_to_normalize]

# Set the font properties
font_properties = fm.FontProperties(family='sans-serif', size=12)

# Plot the heatmap without cell values and remove white borders
plt.figure(figsize=(14, 10))
ax = sns.heatmap(heatmap_weighted_data, annot=False, cmap=custom_cmap, linewidths=0, cbar_kws={'shrink': 0.8})

# Set font properties for labels and title
ax.set_yticklabels(processed_labels, fontproperties=font_properties, fontsize=10)
ax.set_xticklabels(heatmap_weighted_data.columns, rotation=45, ha='right', fontproperties=font_properties, fontsize=10)

# Set axis labels
ax.set_ylabel('binders', fontproperties=font_properties, fontsize=12)

# Set the title
plt.title('Top 50 Binders Based on Weighted Composite Score', fontproperties=font_properties, fontsize=14)

# Adjust layout
plt.tight_layout()

# Save the heatmap as EPS without showing it
plt.savefig('heatmap_top_50.eps', format='eps', bbox_inches='tight')
plt.close()

import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
import matplotlib.font_manager as fm
import argparse

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

# Parse command line arguments
parser = argparse.ArgumentParser(description="Generate heatmap from CSV file.")
parser.add_argument('input_csv', type=str, help="Path to the input CSV file.")
args = parser.parse_args()

# Determine the title based on the input file
if 'minmax_scaler' in args.input_csv:
    title = "Top 50 Binders Based on Weighted Composite Score with Min-Max Scaling Model"
    output_heatmap = 'heatmap_top_50_minmax_scaler.eps'
elif 'standard_scaler' in args.input_csv:
    title = "Top 50 Binders Based on Weighted Composite Score with Standard (Z-score) Scaling Model"
    output_heatmap = 'heatmap_top_50_standard_scaler.eps'
else:
    raise ValueError("Input CSV filename must contain 'minmax_scaler' or 'standard_scaler' to determine the title.")

# Load the processed top 50 data from the provided CSV file
top_50_weighted_data = pd.read_csv(args.input_csv, index_col='description')

# Process y-axis labels to show only the differing text
processed_labels = get_unique_suffixes(top_50_weighted_data.index)

# Define a custom color map similar to the provided image
custom_cmap = sns.color_palette("Blues", as_cmap=True)

# Select the columns used in normalization for the heatmap
metrics_to_normalize = [
    'binder_aligned_rmsd', 
    'pae_binder', 
    'pae_interaction', 
    'plddt_binder', 
    'bsa_score', 
    'salt_bridges', 
    'h_bonds'
]  # Ensure this matches the previous script
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
plt.title(title, fontproperties=font_properties, fontsize=14)

# Adjust layout
plt.tight_layout()

# Save the heatmap as EPS without showing it
plt.savefig(output_heatmap, format='eps', bbox_inches='tight')
plt.close()

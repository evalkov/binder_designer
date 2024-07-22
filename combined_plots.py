
import pandas as pd
import matplotlib.pyplot as plt
import numpy as np

# Load the CSV files
top50_path = 'top50_weighted_composite_scores.csv'
final_path = 'final_binders_list.csv'

top50_df = pd.read_csv(top50_path)
final_df = pd.read_csv(final_path)

# Merging the datasets on the 'description' column
merged_df = pd.merge(top50_df, final_df, on='description', suffixes=('_top50', '_final'))

# Identify the top 50 entries by 'weighted_composite_score'
top_50_indices = merged_df.nlargest(50, 'weighted_composite_score').index

# Function to plot bar plots with specified styles and integer y-axes
def plot_bar_with_integers(column, ax1):
    # Count values for all data
    all_counts = merged_df[column].value_counts().sort_index()
    # Count values for top 50
    top_50_counts = merged_df.loc[top_50_indices, column].value_counts().sort_index()

    # Ensure both Series have the same index
    combined_index = all_counts.index.union(top_50_counts.index)
    all_counts = all_counts.reindex(combined_index, fill_value=0)
    top_50_counts = top_50_counts.reindex(combined_index, fill_value=0)

    # Width for each bar group
    width = 0.35

    # Position of bars on x-axis
    positions = np.arange(len(combined_index))

    # Plot all binders
    ax1.bar(positions - width/2, all_counts, width=width, edgecolor='black', color='#d3d3d3', label='All binders')

    # Create a second y-axis for top 50 binders
    ax2 = ax1.twinx()
    ax2.bar(positions + width/2, top_50_counts, width=width, edgecolor='black', color='red', alpha=0.75, label='Top 50 binders')

    # Set labels and titles
    title = column.replace('_top50', '').replace('_final', '')
    ax1.set_title(f'{title}')
    ax1.set_ylabel('All Binders Frequency')
    ax2.set_ylabel('Top 50 Binders Frequency')

    # Set x-ticks
    ax1.set_xticks(positions)
    ax1.set_xticklabels(combined_index)

    # Set grid lines
    ax1.grid(True, linestyle=':', linewidth=0.5)
    ax2.grid(False)

    # Set y-ticks to integers only
    ax1.yaxis.set_major_locator(plt.MaxNLocator(integer=True))
    ax2.yaxis.set_major_locator(plt.MaxNLocator(integer=True))

    # Combine legends
    lines1, labels1 = ax1.get_legend_handles_labels()
    lines2, labels2 = ax2.get_legend_handles_labels()
    ax1.legend(lines1 + lines2, labels1 + labels2, loc='upper right')

# Exclude 'binder_len', 'salt_bridges', and 'h_bonds' from numeric columns
numeric_columns = merged_df.select_dtypes(include=['float64', 'int64']).columns
numeric_columns_top50 = [col for col in numeric_columns if col.endswith('_top50') and 'binder_len' not in col and 'salt_bridges' not in col and 'h_bonds' not in col]

# Include 'weighted_composite_score' column for scatter plot
numeric_columns_top50.append('weighted_composite_score')

# Number of total plots (2 bar plots + scatter plots for numeric columns)
num_plots = 2 + len(numeric_columns_top50)
fig, axes = plt.subplots((num_plots + 1) // 2, 2, figsize=(20, 4 * ((num_plots + 1) // 2)), sharex=False)
axes = axes.flatten()

# Plot bar plots for 'salt_bridges_final' and 'h_bonds_final' with integer y-axes
plot_bar_with_integers('salt_bridges_final', axes[0])
plot_bar_with_integers('h_bonds_final', axes[1])

# Plotting scatter plots for numeric columns with "_top50" suffix, with top 50 marked in red
for i, col in enumerate(numeric_columns_top50):
    ax = axes[i + 2]
    title = col.replace('_top50', '').replace('_final', '')
    ax.scatter(merged_df.index, merged_df[col], facecolors='none', edgecolors='black', s=10, linewidths=0.5, label='All binders', marker='o')
    ax.scatter(top_50_indices, merged_df.loc[top_50_indices, col], color='red', s=20, label='Top 50 binders', marker='o')
    ax.set_title(f'{title}')
    ax.set_ylabel(f'{title} normalized_score')  # Add y-axis label
    ax.legend()
    ax.grid(True, linestyle=':', linewidth=0.5)

# Turn off any empty subplots
for ax in axes[num_plots:]:
    ax.axis('off')

plt.tight_layout(rect=[0, 0.03, 1, 0.95])

# Save the combined figure as an EPS file
output_combined_plots_path = 'combined_plots.eps'
plt.savefig(output_combined_plots_path, format='eps')
plt.close()

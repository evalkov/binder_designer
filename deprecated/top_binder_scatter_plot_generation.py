import pandas as pd
import matplotlib.pyplot as plt

# Load the dataset
file_path = 'top50_common.csv'
data = pd.read_csv(file_path)

# Plot the two data columns
x_column = 'standard_scale_score'
y_column = 'minmax_scale_score'

# Save the plot as an EPS file
output_path = 'top_binders_scatter_plot.eps'

plt.figure(figsize=(10, 6))
plt.scatter(data[x_column], data[y_column], facecolors='none', edgecolors='black', s=50, linewidths=0.5, marker='o')
plt.title('Scatter Plot of Top-Scoring Binders')
plt.xlabel('Standard (Z-score) Scale Score')
plt.ylabel('Min-Max Scale Score')
plt.grid(True, linestyle=':', linewidth=0.5)

plt.savefig(output_path, format='eps')
plt.close()

print(f"The scatter plot has been saved as {output_path}.")

import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns

# Load the final_binders_list.csv
final_binders_list = pd.read_csv('final_binders_list.csv')

# Compute the correlation matrix
corr_matrix = final_binders_list.corr(numeric_only=True)

# Display the correlation matrix
plt.figure(figsize=(12, 10))
sns.heatmap(corr_matrix, annot=True, cmap='coolwarm', linewidths=0.5)
plt.title('Correlation Matrix')

# Adjust the plot to avoid truncated labels
plt.xticks(rotation=45, ha='right')
plt.yticks(rotation=0)
plt.tight_layout()

# Save the plot as an EPS file
plt.savefig('correlation_matrix.eps', format='eps')

# Show the plot
#plt.show()

import pandas as pd
import matplotlib.pyplot as plt
from sklearn.cluster import KMeans
from sklearn.metrics import silhouette_score
import os

def analyze_binders(csv_file_path):
    # Load the CSV file
    binders_df = pd.read_csv(csv_file_path)
    
    # Selecting the relevant features for clustering
    features = binders_df[['bsa_score', 'salt_bridges']].copy()
    
    # Drop rows with missing values in 'bsa_score'
    features_cleaned = features.dropna().copy()
    
    # Perform k-means clustering on the cleaned data
    kmeans = KMeans(n_clusters=3, random_state=42)
    features_cleaned.loc[:, 'cluster'] = kmeans.fit_predict(features_cleaned)
    
    # Calculate the silhouette score for the clustering
    silhouette_avg = silhouette_score(features_cleaned, features_cleaned['cluster'])
    print(f"Silhouette Score: {silhouette_avg}")
    
    # Adding cluster labels back to the main dataframe
    binders_df_cleaned = binders_df.dropna(subset=['bsa_score']).copy()
    binders_df_cleaned.loc[features_cleaned.index, 'cluster'] = features_cleaned['cluster']
    
    # Summary of each cluster
    cluster_summary = binders_df_cleaned.groupby('cluster').agg({
        'bsa_score': ['mean', 'std', 'min', 'max', 'count'],
        'salt_bridges': ['mean', 'std', 'min', 'max'],
        'binder_aligned_rmsd': ['mean', 'std'],
        'pae_binder': ['mean', 'std'],
        'pae_interaction': ['mean', 'std'],
        'plddt_binder': ['mean', 'std'],
        'h_bonds': ['mean', 'std'],
    })
    
    print("Cluster Summary:")
    print(cluster_summary)
    
    # Identifying the top 20 scoring binders based on bsa_score
    top_20_binders = binders_df_cleaned.nlargest(20, 'bsa_score')
    
    # Writing the top 20 binders' stats to a text file
    with open('top20.txt', 'w') as f:
        f.write("Top 20 Scoring Binders:\n")
        for index, row in top_20_binders.iterrows():
            f.write(f"{row['description']}\n")
            f.write(f"BSA Score: {row['bsa_score']}\n")
            f.write(f"Salt Bridges: {row['salt_bridges']}\n")
            f.write(f"H Bonds: {row['h_bonds']}\n")
            f.write(f"Binder Aligned RMSD: {row['binder_aligned_rmsd']}\n")
            f.write(f"PAE Binder: {row['pae_binder']}\n")
            f.write(f"PAE Interaction: {row['pae_interaction']}\n")
            f.write(f"PLDDT Binder: {row['plddt_binder']}\n\n")
    
    return binders_df_cleaned, top_20_binders

def plot_salt_bridges_vs_bsa(binders_df_cleaned, top_20_binders):
    # Scatter plot of Salt Bridges vs BSA Score
    plt.figure(figsize=(10, 6))
    plt.scatter(binders_df_cleaned['bsa_score'], binders_df_cleaned['salt_bridges'], c=binders_df_cleaned['cluster'], cmap='viridis', label='Cluster')
    
    # Marking the top 20 binders
    plt.scatter(top_20_binders['bsa_score'], top_20_binders['salt_bridges'], c='red', edgecolor='k', s=150, label='Top 20 Binders')
    
    plt.colorbar(label='Cluster')
    plt.xlabel('BSA Score')
    plt.ylabel('Salt Bridges')
    plt.title('Salt Bridges vs BSA Score')
    plt.legend()
    plt.savefig('plot_salt_bridges_vs_bsa.png', format='png')
    plt.show()

def plot_pae_vs_bsa(binders_df_cleaned, top_20_binders):
    # Scatter plot of PAE Interaction vs BSA Score
    plt.figure(figsize=(10, 6))
    plt.scatter(binders_df_cleaned['bsa_score'], binders_df_cleaned['pae_interaction'], c=binders_df_cleaned['cluster'], cmap='viridis', label='Cluster')
    
    # Marking the top 20 binders
    plt.scatter(top_20_binders['bsa_score'], top_20_binders['pae_interaction'], c='red', edgecolor='k', s=150, label='Top 20 Binders')
    
    plt.colorbar(label='Cluster')
    plt.xlabel('BSA Score')
    plt.ylabel('PAE Interaction')
    plt.title('PAE Interaction vs BSA Score')
    plt.legend()
    plt.savefig('plot_pae_vs_bsa.png', format='png')
    plt.show()

def plot_plddt_vs_bsa(binders_df_cleaned, top_20_binders):
    # Scatter plot of PLDDT Binder vs BSA Score
    plt.figure(figsize=(10, 6))
    plt.scatter(binders_df_cleaned['bsa_score'], binders_df_cleaned['plddt_binder'], c=binders_df_cleaned['cluster'], cmap='viridis', label='Cluster')
    
    # Marking the top 20 binders
    plt.scatter(top_20_binders['bsa_score'], top_20_binders['plddt_binder'], c='red', edgecolor='k', s=150, label='Top 20 Binders')
    
    plt.colorbar(label='Cluster')
    plt.xlabel('BSA Score')
    plt.ylabel('PLDDT Binder')
    plt.title('PLDDT Binder vs BSA Score')
    plt.legend()
    plt.savefig('plot_plddt_vs_bsa.png', format='png')
    plt.show()

def plot_rmsd_vs_bsa(binders_df_cleaned, top_20_binders):
    # Scatter plot of Binder Aligned RMSD vs BSA Score
    plt.figure(figsize=(10, 6))
    plt.scatter(binders_df_cleaned['bsa_score'], binders_df_cleaned['binder_aligned_rmsd'], c=binders_df_cleaned['cluster'], cmap='viridis', label='Cluster')
    
    # Marking the top 20 binders
    plt.scatter(top_20_binders['bsa_score'], top_20_binders['binder_aligned_rmsd'], c='red', edgecolor='k', s=150, label='Top 20 Binders')
    
    plt.colorbar(label='Cluster')
    plt.xlabel('BSA Score')
    plt.ylabel('Binder Aligned RMSD')
    plt.title('Binder Aligned RMSD vs BSA Score')
    plt.legend()
    plt.savefig('plot_rmsd_vs_bsa.png', format='png')
    plt.show()

def generate_cxc_file(csv_file_path, output_cxc_path):
    binders_df_cleaned, top_20_binders = analyze_binders(csv_file_path)
    
    # Generating the binders.cxc content
    cxc_content = "set bgcolor white\n"
    for index, binder in enumerate(top_20_binders['description']):
        cxc_content += f"open {binder}.pdb\n"
    cxc_content += (
        "cartoon style protein modeh tube rad 2 sides 24\n"
        "cartoon style width 2 thick 0.2\n"
        "rainbow chain palette RdYlBu-5\n"
        "lighting simple shadows false intensity 0.5\n"
        "view all\n"
        "hide atoms\n"
        "show cartoons\n"
        "hide all models\n"
        "show #1 models\n"
        "matchmaker all to #1/B pairing bs\n"
    )
    
    for i in range(1, 21):
        cxc_content += (
            f"interfaces select #{i}/B contacting #{i}/A bothSides true\n"
            f"contacts #{i}/A restrict #{i}/B intraMol false\n"
            "show sel atoms\n"
            "select clear\n"
        )
    
    cxc_content += (
        "delete H\n"
        "color byhetero\n"
    )
    
    # Ensure the output directory exists
    os.makedirs(os.path.dirname(output_cxc_path), exist_ok=True)
    
    # Writing the cxc content to the output file
    with open(output_cxc_path, 'w') as file:
        file.write(cxc_content)
    
    print(f"Generated binders.cxc script and saved to {output_cxc_path}")
    
    # Plotting clusters and top binders
    plot_salt_bridges_vs_bsa(binders_df_cleaned, top_20_binders)
    plot_pae_vs_bsa(binders_df_cleaned, top_20_binders)
    plot_plddt_vs_bsa(binders_df_cleaned, top_20_binders)
    plot_rmsd_vs_bsa(binders_df_cleaned, top_20_binders)

# Example usage
csv_file_path = '/Users/valkove2/Downloads/Archive/final_binders_list.csv'  # Path to the CSV file
output_cxc_path = '/Users/valkove2/Downloads/Archive/top10/binders.cxc'  # Path where the generated binders.cxc file will be saved
generate_cxc_file(csv_file_path, output_cxc_path)

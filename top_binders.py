import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
import re  # Import regular expressions module

def merge_csv_files(standard_scaler_path, minmax_scaler_path, output_csv_path):
    """
    Merges two CSV files on the 'description' column and saves the merged DataFrame to a new CSV file.

    Parameters:
    - standard_scaler_path (str): Path to the standard scaler CSV file.
    - minmax_scaler_path (str): Path to the min-max scaler CSV file.
    - output_csv_path (str): Path where the merged CSV will be saved.
    """
    # Load the CSV files
    try:
        df1 = pd.read_csv(standard_scaler_path)
        df2 = pd.read_csv(minmax_scaler_path)
        print("CSV files loaded successfully.")
    except FileNotFoundError as e:
        print(f"Error: {e}")
        return
    except pd.errors.EmptyDataError:
        print("Error: One of the CSV files is empty.")
        return
    except pd.errors.ParserError:
        print("Error: Parsing error while reading the CSV files.")
        return

    # Print column names for verification
    print("\nColumns in df1:", df1.columns.tolist())
    print("Columns in df2:", df2.columns.tolist())

    # Standardize column names by stripping spaces and converting to lowercase
    df1.columns = df1.columns.str.strip().str.lower()
    df2.columns = df2.columns.str.strip().str.lower()

    # Verify that 'description' exists in both DataFrames
    if 'description' not in df1.columns:
        raise KeyError("Column 'description' not found in the standard scaler CSV.")
    if 'description' not in df2.columns:
        raise KeyError("Column 'description' not found in the min-max scaler CSV.")

    # Check for duplicates in 'description'
    duplicates_df1 = df1[df1.duplicated(subset=['description'], keep=False)]
    duplicates_df2 = df2[df2.duplicated(subset=['description'], keep=False)]

    if not duplicates_df1.empty:
        print("\nWarning: Duplicates found in df1 'description' column:")
        print(duplicates_df1[['description']])
        # Handle duplicates as needed

    if not duplicates_df2.empty:
        print("\nWarning: Duplicates found in df2 'description' column:")
        print(duplicates_df2[['description']])
        # Handle duplicates as needed

    # Merge the dataframes on the 'description' column
    merged_df = pd.merge(df1, df2, on='description', suffixes=('_standard', '_minmax'))
    print("\nDataFrames merged successfully.")

    # Select relevant columns and rename them
    try:
        merged_df = merged_df[['description',
                               'sequence_standard',  # Assuming 'sequence' exists in df1
                               'weighted_composite_score_standard',
                               'weighted_composite_score_minmax']]
    except KeyError as e:
        print(f"Error: {e}")
        print("Available columns after merge:", merged_df.columns.tolist())
        return

    merged_df.columns = ['description', 'sequence', 'standard_scale_score', 'minmax_scale_score']

    # Handle missing values if any in the scores
    if merged_df[['standard_scale_score', 'minmax_scale_score']].isnull().any().any():
        print("\nWarning: Missing values detected in the score columns. Dropping rows with missing values.")
        merged_df = merged_df.dropna(subset=['standard_scale_score', 'minmax_scale_score'])

    # Save the merged dataframe to a new CSV file
    merged_df.to_csv(output_csv_path, index=False)
    print(f"\nMerged CSV saved successfully as '{output_csv_path}'.")

    return merged_df

def extract_label(description):
    """
    Extracts the first two groups of digits from a description string and joins them with an underscore.

    Parameters:
    - description (str): The description string to process.

    Returns:
    - label (str): The extracted label in the format 'number1_number2'.
    """
    numbers = re.findall(r'\d+', description)
    if len(numbers) >= 2:
        label = numbers[0] + '_' + numbers[1]
    elif len(numbers) == 1:
        label = numbers[0]
    else:
        label = ''
    return label

def generate_scatter_plot(data, x_column, y_column, output_plot_path):
    """
    Generates a scatter plot for two specified columns, labels the top 5 data points with adjusted positions
    to avoid overlap, and saves it as an EPS file.

    Parameters:
    - data (pd.DataFrame): DataFrame containing the data to plot.
    - x_column (str): Column name for the x-axis.
    - y_column (str): Column name for the y-axis.
    - output_plot_path (str): Path where the plot will be saved.
    """
    import numpy as np

    # Verify that the specified columns exist in the DataFrame
    if x_column not in data.columns:
        raise KeyError(f"Column '{x_column}' not found in the DataFrame.")
    if y_column not in data.columns:
        raise KeyError(f"Column '{y_column}' not found in the DataFrame.")

    # Extract labels from the 'description' column
    data['label'] = data['description'].apply(extract_label)

    # Identify the top 5 data points based on 'standard_scale_score'
    top_n = 5
    top_indices = data['standard_scale_score'].nlargest(top_n).index

    # Prepare colors: light grey for all points, red for top 5
    colors = ['lightgrey'] * len(data)
    for idx in top_indices:
        colors[idx] = 'red'

    # Plot the data points
    plt.figure(figsize=(12, 8))
    plt.scatter(data[x_column], data[y_column],
                facecolors=colors, edgecolors='black',
                s=50, linewidths=0.5, marker='o')

    plt.title('Scatter Plot of Top-Scoring Binders')
    plt.xlabel('Standard (Z-score) Scale Score')
    plt.ylabel('Min-Max Scale Score')
    plt.grid(True, linestyle=':', linewidth=0.5)

    # Add labels to the top 5 points with adjusted positions
    x = data[x_column].values
    y = data[y_column].values
    labels = data['label'].values

    # Adjust label font size
    label_font_size = 12  # Larger labels

    # Adjust offset to bring labels a bit further from data points
    offset = 0.02  # Increased offset to move labels slightly further away
    xytext_factor = 600  # Adjusted to control label line length

    # Initialize arrays for offsets
    x_offsets = np.zeros(len(x))
    y_offsets = np.zeros(len(y))

    # Keep track of used positions to minimize overlap
    used_positions = []

    for idx in top_indices:
        i = idx  # Current index

        # Try different positions around the data point
        positions = [
            (offset, offset),
            (-offset, offset),
            (-offset, -offset),
            (offset, -offset),
            (0, offset),
            (0, -offset),
            (offset, 0),
            (-offset, 0),
        ]

        for x_off, y_off in positions:
            proposed_pos = (x[i] + x_off, y[i] + y_off)
            overlap = False
            for pos in used_positions:
                if abs(proposed_pos[0] - pos[0]) < offset/2 and abs(proposed_pos[1] - pos[1]) < offset/2:
                    overlap = True
                    break
            if not overlap:
                x_offsets[i] = x_off
                y_offsets[i] = y_off
                used_positions.append(proposed_pos)
                break
        else:
            # Default position if all positions overlap
            x_offsets[i] = offset
            y_offsets[i] = offset
            used_positions.append((x[i] + offset, y[i] + offset))

        # Adjust arrow properties to prevent penetration into data point
        arrowprops = dict(
            arrowstyle='-',
            color='gray',
            lw=0.5,
            shrinkA=0,
            shrinkB=5
        )

        plt.annotate(
            labels[i],
            (x[i], y[i]),
            textcoords="offset points",
            xytext=(x_offsets[i]*xytext_factor, y_offsets[i]*xytext_factor),
            ha='center',
            fontsize=label_font_size,
            arrowprops=arrowprops
        )

    # Adjust layout for better spacing
    plt.tight_layout()

    # Save the plot as an EPS file
    plt.savefig(output_plot_path, format='eps')
    plt.close()

    print(f"The scatter plot has been saved as '{output_plot_path}'.")

def main():
    # Define file paths
    standard_scaler_csv = 'top_50_binders_weighted_standard_scaler.csv'
    minmax_scaler_csv = 'top_50_binders_weighted_minmax_scaler.csv'
    merged_csv = 'top50_common.csv'
    scatter_plot_eps = 'top_binders_scatter_plot.eps'

    # Step 1: Merge the CSV files
    merged_data = merge_csv_files(standard_scaler_csv, minmax_scaler_csv, merged_csv)

    if merged_data is None:
        print("Merging CSV files failed. Exiting script.")
        return

    # Optional: Display the first few rows of the merged DataFrame
    print("\nFirst few rows of the merged DataFrame:")
    print(merged_data.head())

    # Step 2: Generate the scatter plot with adjusted labels
    generate_scatter_plot(merged_data,
                          x_column='standard_scale_score',
                          y_column='minmax_scale_score',
                          output_plot_path=scatter_plot_eps)

if __name__ == "__main__":
    main()


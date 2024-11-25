#!/usr/bin/env python3

import pandas as pd
import os
import argparse
import logging
from pathlib import Path
import sys

def setup_logging(log_file: Path):
    """
    Set up logging to output to both console and a log file.

    Args:
        log_file (Path): Path to the log file.
    """
    logging.basicConfig(
        level=logging.INFO,
        format='%(asctime)s [%(levelname)s] %(message)s',
        handlers=[
            logging.FileHandler(log_file),
            logging.StreamHandler(sys.stdout)
        ]
    )

def parse_arguments():
    """
    Parse command-line arguments.

    Returns:
        argparse.Namespace: Parsed arguments.
    """
    parser = argparse.ArgumentParser(description="Merge binders.csv with contacts.csv to create final_binders_list.csv.")
    parser.add_argument('--binders_csv', type=str, default='binders.csv', help='Path to the binders.csv file.')
    parser.add_argument('--contacts_csv', type=str, default='contacts.csv', help='Path to the contacts.csv file.')
    parser.add_argument('--output_csv', type=str, default='final_binders_list.csv', help='Path to the output CSV file.')
    parser.add_argument('--log_file', type=str, default='merge_binders_contacts.log', help='Path to the log file.')
    return parser.parse_args()

def check_file_exists(file_path: Path):
    """
    Check if a file exists. Exit the script if it does not.

    Args:
        file_path (Path): Path to the file.
    """
    if not file_path.is_file():
        logging.error(f"File not found: {file_path}")
        sys.exit(1)
    else:
        logging.info(f"Found file: {file_path}")

def read_csv_file(file_path: Path, file_description: str) -> pd.DataFrame:
    """
    Read a CSV file into a pandas DataFrame.

    Args:
        file_path (Path): Path to the CSV file.
        file_description (str): Description of the file for logging.

    Returns:
        pd.DataFrame: The loaded DataFrame.
    """
    try:
        df = pd.read_csv(file_path)
        logging.info(f"Successfully read {file_description} from {file_path}")
        return df
    except Exception as e:
        logging.error(f"Failed to read {file_description} from {file_path}: {e}")
        sys.exit(1)

def validate_columns(df: pd.DataFrame, required_columns: list, df_description: str):
    """
    Validate that the required columns exist in the DataFrame.

    Args:
        df (pd.DataFrame): The DataFrame to validate.
        required_columns (list): List of required column names.
        df_description (str): Description of the DataFrame for logging.
    """
    missing_columns = [col for col in required_columns if col not in df.columns]
    if missing_columns:
        logging.error(f"Missing columns in {df_description}: {missing_columns}")
        sys.exit(1)
    else:
        logging.info(f"All required columns found in {df_description}")

def merge_dataframes(binders_df: pd.DataFrame, contacts_df: pd.DataFrame) -> pd.DataFrame:
    """
    Merge binders_df with contacts_df on binder identifiers.

    Args:
        binders_df (pd.DataFrame): DataFrame containing binder information.
        contacts_df (pd.DataFrame): DataFrame containing contact information.

    Returns:
        pd.DataFrame: The merged DataFrame.
    """
    # Assuming 'description' in binders_df corresponds to 'binder' in contacts_df
    try:
        combined_df = pd.merge(
            binders_df,
            contacts_df[['binder', 'bsa_score', 'salt_bridges', 'h_bonds']],
            left_on='description',
            right_on='binder',
            how='inner'
        )
        logging.info("Successfully merged binders.csv with contacts.csv")
        return combined_df
    except KeyError as e:
        logging.error(f"Merge failed due to missing key: {e}")
        sys.exit(1)
    except Exception as e:
        logging.error(f"An unexpected error occurred during merging: {e}")
        sys.exit(1)

def clean_combined_dataframe(combined_df: pd.DataFrame) -> pd.DataFrame:
    """
    Clean the merged DataFrame by dropping duplicate columns.

    Args:
        combined_df (pd.DataFrame): The merged DataFrame.

    Returns:
        pd.DataFrame: The cleaned DataFrame.
    """
    if 'binder' in combined_df.columns:
        combined_df.drop(columns=['binder'], inplace=True)
        logging.info("Dropped duplicate 'binder' column after merging")
    else:
        logging.warning("'binder' column not found in merged DataFrame. Skipping drop.")
    return combined_df

def save_to_csv(df: pd.DataFrame, output_path: Path):
    """
    Save the DataFrame to a CSV file.

    Args:
        df (pd.DataFrame): The DataFrame to save.
        output_path (Path): Path to the output CSV file.
    """
    try:
        df.to_csv(output_path, index=False)
        logging.info(f"Successfully saved merged data to {output_path}")
    except Exception as e:
        logging.error(f"Failed to save merged data to {output_path}: {e}")
        sys.exit(1)

def main():
    # Parse command-line arguments
    args = parse_arguments()
    binders_csv_path = Path(args.binders_csv)
    contacts_csv_path = Path(args.contacts_csv)
    output_csv_path = Path(args.output_csv)
    log_file_path = Path(args.log_file)

    # Set up logging
    setup_logging(log_file_path)

    logging.info("=== Starting binder and contacts merge process ===")

    # Check if input files exist
    check_file_exists(binders_csv_path)
    check_file_exists(contacts_csv_path)

    # Read input CSV files
    binders_df = read_csv_file(binders_csv_path, "binders.csv")
    contacts_df = read_csv_file(contacts_csv_path, "contacts.csv")

    # Validate required columns in binders.csv (all lowercase)
    required_binders_columns = ["binder_aligned_rmsd", "pae_binder", "pae_interaction", "plddt_binder", "description", "sequence", "size"]
    validate_columns(binders_df, required_binders_columns, "binders.csv")

    # Validate required columns in contacts.csv
    required_contacts_columns = ["binder", "bsa_score", "salt_bridges", "h_bonds"]
    validate_columns(contacts_df, required_contacts_columns, "contacts.csv")

    # Merge the DataFrames
    combined_df = merge_dataframes(binders_df, contacts_df)

    # Clean the merged DataFrame
    combined_df = clean_combined_dataframe(combined_df)

    # Save the merged DataFrame to final_binders_list.csv
    save_to_csv(combined_df, output_csv_path)

    logging.info("=== Binder and contacts merge process completed successfully ===")

if __name__ == "__main__":
    main()


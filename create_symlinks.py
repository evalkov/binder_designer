#!/usr/bin/env python3

import os
import argparse
import logging
from pathlib import Path
import sys
import concurrent.futures
import threading

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
    parser = argparse.ArgumentParser(
        description="Create symbolic links for .pdb files from a source directory to the target directory."
    )
    parser.add_argument(
        '--source_dir', '-s',
        type=str,
        required=True,
        help='Path to the source directory containing .pdb files.'
    )
    parser.add_argument(
        '--target_dir', '-t',
        type=str,
        default=os.getcwd(),
        help='Path to the target directory where symlinks will be created. Defaults to the current working directory.'
    )
    parser.add_argument(
        '--settings_file', '-c',
        type=str,
        required=True,
        help='Path to the settings file containing configuration parameters.'
    )
    parser.add_argument(
        '--log_file', '-l',
        type=str,
        default='create_symlinks.log',
        help='Path to the log file.'
    )
    return parser.parse_args()

def parse_settings_file(settings_file: Path) -> dict:
    """
    Parse the settings file to extract configuration parameters.

    Args:
        settings_file (Path): Path to the settings file.

    Returns:
        dict: Dictionary containing configuration parameters.
    """
    if not settings_file.is_file():
        logging.error(f"Settings file does not exist: {settings_file}")
        sys.exit(1)
    
    settings = {}
    try:
        with open(settings_file, 'r') as f:
            for line in f:
                # Ignore empty lines and comments
                line = line.strip()
                if not line or line.startswith('#'):
                    continue
                if ':' in line:
                    key, value = line.split(':', 1)
                    settings[key.strip()] = value.strip()
    except Exception as e:
        logging.error(f"Failed to read settings file {settings_file}: {e}")
        sys.exit(1)
    
    # Validate required settings
    if 'max_workers' not in settings:
        logging.error(f"'max_workers' not found in settings file {settings_file}.")
        sys.exit(1)
    
    try:
        settings['max_workers'] = int(settings['max_workers'])
    except ValueError:
        logging.error(f"'max_workers' must be an integer in settings file {settings_file}.")
        sys.exit(1)
    
    # Optional settings can be handled here
    # For example:
    # settings['batch_size'] = int(settings.get('batch_size', 100))
    
    logging.info(f"Settings loaded from {settings_file}: {settings}")
    return settings

def validate_directories(source_dir: Path, target_dir: Path):
    """
    Validate that source and target directories exist.

    Args:
        source_dir (Path): Source directory path.
        target_dir (Path): Target directory path.
    """
    if not source_dir.is_dir():
        logging.error(f"Source directory does not exist or is not a directory: {source_dir}")
        sys.exit(1)
    else:
        logging.info(f"Source directory validated: {source_dir}")

    if not target_dir.exists():
        try:
            target_dir.mkdir(parents=True, exist_ok=True)
            logging.info(f"Target directory created: {target_dir}")
        except Exception as e:
            logging.error(f"Failed to create target directory {target_dir}: {e}")
            sys.exit(1)
    elif not target_dir.is_dir():
        logging.error(f"Target path exists but is not a directory: {target_dir}")
        sys.exit(1)
    else:
        logging.info(f"Target directory validated: {target_dir}")

def create_symlink(file_path: Path, target_dir: Path, lock: threading.Lock, counters: dict):
    """
    Create a symbolic link for the given file in the target directory.

    Args:
        file_path (Path): Path to the source file.
        target_dir (Path): Path to the target directory.
        lock (threading.Lock): A lock object to synchronize access to counters.
        counters (dict): A dictionary to keep track of created symlinks.
    """
    target_path = target_dir / file_path.name
    try:
        if not target_path.exists():
            os.symlink(file_path.resolve(), target_path)
            with lock:
                counters['created_symlinks'] += 1
            logging.info(f"Symlink created: {target_path} -> {file_path.resolve()}")
        else:
            logging.warning(f"Symlink already exists: {target_path}")
    except Exception as e:
        logging.error(f"Failed to create symlink for {file_path}: {e}")

def main():
    # Parse command-line arguments
    args = parse_arguments()
    source_dir = Path(args.source_dir)
    target_dir = Path(args.target_dir)
    settings_file = Path(args.settings_file)
    log_file = Path(args.log_file)

    # Set up logging
    setup_logging(log_file)
    logging.info("=== Starting symlink creation process ===")

    # Parse settings file
    settings = parse_settings_file(settings_file)
    max_workers = settings['max_workers']

    # Validate directories
    validate_directories(source_dir, target_dir)

    # Get a list of all .pdb files in the source directory
    pdb_files = list(source_dir.glob("*.pdb"))
    total_files = len(pdb_files)

    if not pdb_files:
        logging.warning(f"No .pdb files found in the source directory '{source_dir}'. Exiting.")
        sys.exit(0)

    logging.info(f"Found {total_files} .pdb files in '{source_dir}'.")
    logging.info(f"Creating symlinks in '{target_dir}' using up to {max_workers} workers...")

    # Initialize counters and lock for thread safety
    counters = {'created_symlinks': 0}
    lock = threading.Lock()

    # Use ThreadPoolExecutor for parallel symlink creation
    with concurrent.futures.ThreadPoolExecutor(max_workers=max_workers) as executor:
        futures = [
            executor.submit(create_symlink, file_path, target_dir, lock, counters)
            for file_path in pdb_files
        ]
        # Optionally, monitor progress
        for idx, future in enumerate(concurrent.futures.as_completed(futures), 1):
            try:
                future.result()
            except Exception as e:
                logging.error(f"Unhandled exception: {e}")
            if idx % 1000 == 0 or idx == total_files:
                logging.info(f"Progress: {idx}/{total_files} symlinks processed.")

    # Display summary
    logging.info(f"\nFinished! {counters['created_symlinks']} symlinks were created out of {total_files} files.")
    logging.info("=== Symlink creation process completed ===")

if __name__ == "__main__":
        main()


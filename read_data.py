import os
import pandas as pd
import numpy as np
from tqdm import tqdm
import re
from pathlib import Path
import matplotlib.pyplot as plt
import seaborn as sns

"""
Read all FISH quantification data from a folder and store in a DataFrame.
This program:
1. Recursively finds all FISH quant files in the specified folder
2. Reads each file and extracts relevant data
3. Combines all data into a single DataFrame with metadata
4. Generates summary visualizations
"""

# Configuration
FISH_QUANT_DIR = os.path.join(os.path.dirname(__file__), "fish_quant")
OUTPUT_DIR = os.path.join(os.path.dirname(__file__), "fish_quant_analysis")
FILE_EXTENSIONS = [".csv", ".txt"]

# Create output directory if it doesn't exist
os.makedirs(OUTPUT_DIR, exist_ok=True)

def find_fish_quant_files(dir_path):
    """Find all FISH quantification files in the directory."""
    all_files = []
    
    print(f"Searching for FISH quant files in: {dir_path}")
    
    for root, _, files in os.walk(dir_path):
        for file in files:
            if any(file.lower().endswith(ext) for ext in FILE_EXTENSIONS):
                all_files.append(os.path.join(root, file))
    
    print(f"Found {len(all_files)} potential FISH quant files.")
    return all_files

def extract_metadata_from_path(file_path):
    """Extract metadata from file path and name."""
    filename = os.path.basename(file_path)
    
    # Default metadata
    metadata = {
        "file_name": filename,
        "full_path": file_path
    }
    
    # Extract time points from filename (e.g., -4h-, -8h-)
    time_match = re.search(r"-(\d+)h-", filename)
    if time_match:
        metadata["time_point"] = int(time_match.group(1))
    
    # Check for reinduction condition
    if "Reinduction" in filename:
        metadata["condition"] = "reinduction"
    else:
        metadata["condition"] = "primary"
    
    return metadata

def read_fish_quant_file(file_path):
    """Read and parse FISH quantification file."""
    metadata = extract_metadata_from_path(file_path)
    
    try:
        # Try different delimiters
        for delimiter in [',', '\t', ';', ' ']:
            try:
                data = pd.read_csv(file_path, delimiter=delimiter)
                if not data.empty and len(data.columns) >= 2:
                    metadata["format"] = "delimited"
                    metadata["delimiter"] = delimiter
                    return data, metadata
            except:
                continue
        
        # Fallback to automatic delimiter detection
        data = pd.read_csv(file_path, sep=None, engine='python')
        metadata["format"] = "auto-detected"
        return data, metadata
        
    except Exception as e:
        print(f"Error reading file {file_path}: {e}")
        return None, metadata

def process_fish_quant_data(data, metadata):
    """Process the FISH quantification data and add metadata."""
    if data is None:
        return pd.DataFrame()
    
    # Initialize processed DataFrame
    processed = data.copy()
    
    # Add time point and condition columns
    processed['time_point'] = metadata.get('time_point')
    processed['condition'] = metadata.get('condition', 'unknown')
    
    # Rename data columns based on their content
    rename_dict = {}
    for col in processed.columns:
        if col not in ['time_point', 'condition']:
            col_lower = col.lower()
            if any(term in col_lower for term in ['count', 'spot', 'intens', 'signal']):
                rename_dict[col] = f"value_{col}"
    
    processed.rename(columns=rename_dict, inplace=True)
    return processed

def main():
    """Main function to run the data processing pipeline."""
    # Find all FISH quant files
    quant_files = find_fish_quant_files(FISH_QUANT_DIR)
    
    if not quant_files:
        print(f"No FISH quant files found in {FISH_QUANT_DIR}")
        return pd.DataFrame()
    
    # Process files with progress bar
    all_data = pd.DataFrame()
    
    for file_path in tqdm(quant_files, desc="Reading files"):
        data, metadata = read_fish_quant_file(file_path)
        processed_data = process_fish_quant_data(data, metadata)
        
        if not processed_data.empty:
            all_data = pd.concat([all_data, processed_data], ignore_index=True)
    
    # Save results grouped by time point
    if not all_data.empty:
        time_points = sorted(all_data['time_point'].dropna().unique())
        conditions = all_data['condition'].unique()
        
        # Save data for each time point separately
        for tp in time_points:
            for cond in conditions:
                subset = all_data[
                    (all_data['time_point'] == tp) & 
                    (all_data['condition'] == cond)
                ]
                
                # Group by cell and region, sum values to get counts
                grouped = subset.groupby(['cell', 'frame', 'type']).size().reset_index(name='count')
                
                # Save the count data to CSV file
                output_path = os.path.join(OUTPUT_DIR, f"counts_{int(tp)}h_{cond}.csv")
                grouped.to_csv(output_path, index=False)
                
                # Plot histogram of counts for each value column (using Transcription Site data)
                value_columns = grouped
                
                if not value_columns.empty and not subset.empty:
                    # Create histogram
                    print(f"Creating histogram for {int(tp)}h ({cond})")
                    fig, ax = plt.subplots(figsize=(10, 6))
                    
                    sns.histplot(value_columns, 
                                 binwidth=1, 
                                 alpha=0.7,
                                 stat="density",
                                 label=f"Counts at {int(tp)}h ({cond})")
                    
                    ax.set_title(f"Histogram of FISH counts at {int(tp)}h ({cond})")
                    ax.set_xlabel("Count value")
                    ax.set_ylabel("Frequency")
                    ax.legend()
                    ax.set_xlim(0,200)
                    
                    # Save the histogram
                    plt.tight_layout()
                    plt.savefig(os.path.join(OUTPUT_DIR, f"histogram/{cond}_{int(tp)}h.png"))
                    plt.close()
        
        print("\nData grouped by time points:")
        for tp in time_points:
            n_primary = len(all_data[
                (all_data['time_point'] == tp) & 
                (all_data['condition'] == 'primary')
            ])
            n_reind = len(all_data[
                (all_data['time_point'] == tp) & 
                (all_data['condition'] == 'reinduction')
            ])
            print(f"- {int(tp)}h: Primary={n_primary}, Reinduction={n_reind} samples")
    
    return all_data

if __name__ == "__main__":
    all_data = main()
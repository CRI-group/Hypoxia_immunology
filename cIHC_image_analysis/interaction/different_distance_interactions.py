"""
====================================================================================================
Modified Cell Interaction Calculation for All Pairwise Cell Type Combinations - Grouped by Hypoxia Class
====================================================================================================
Summary:
This script performs the following steps:
  1. Loads the DataFrame with a structure where cell types/markers are in a single column
  2. Reduces the DataFrame to the necessary columns including spatial coordinates, spot identifiers,
     and the cell type/marker column
  3. For a specified maximum distance, calculates interaction indices for ALL PAIRWISE COMBINATIONS
     of cell types, GROUPED BY BOTH IMAGE ID AND HYPOXIA CLASS
  4. The interaction index is computed per image-hypoxia class combination as:
         (Number of reference cells with at least one nearby target cell)
         divided by the square root of (total reference cell count × total target cell count)
  5. Saves the final aggregated DataFrame with all interaction indices to a CSV file
====================================================================================================
"""

import os
import sys
import numpy as np
import pandas as pd
from scipy.spatial import cKDTree
import logging
import warnings

# Configure logging
logger = logging.getLogger()
logger.setLevel(logging.INFO)
if not logger.handlers:
    handler = logging.StreamHandler(sys.stdout)
    formatter = logging.Formatter('%(asctime)s - %(levelname)s - %(message)s')
    handler.setFormatter(formatter)
    logger.addHandler(handler)
else:
    # Add a handler if not already added
    handler = logging.StreamHandler(sys.stdout)
    handler.setFormatter(logging.Formatter('%(asctime)s - %(levelname)s - %(message)s'))
    logger.addHandler(handler)

logging.info("Starting modified cell interaction calculation with all pairwise combinations...")

# ----------------------------------------------------------------------
# Step 1: Load data or verify it's loaded
# ----------------------------------------------------------------------
try:
    # Uncomment and adjust the following lines if you need to load the DataFrame
    input_df_path = "/scratch/svc_td_cri/projects/multiplex/Old_narvi/hypoxia/intensities/combined_Hypoxia_alldata.csv"
    df = pd.read_csv(input_df_path, low_memory=False)
    logging.info(f"DataFrame loaded from {input_df_path}. Shape: {df.shape}")
except Exception as e:
    logging.error(f"Error loading DataFrame: {str(e)}")
    raise

# ----------------------------------------------------------------------
# Step 2: Identify necessary columns
# ----------------------------------------------------------------------
columns_to_keep = [
    'Identifier',
    'Centroid.X.um',
    'Centroid.Y.um',
    'Parent',
    'Classification'  # This is the column containing your cell type/marker values
]

# Check if all necessary columns exist in the dataframe
missing_columns = [col for col in columns_to_keep if col not in df.columns]
if missing_columns:
    logging.error(f"Missing columns in dataframe: {missing_columns}")
    raise ValueError(f"Missing columns in dataframe: {missing_columns}")

reduced_df = df[columns_to_keep].copy()
logging.info(f"Columns reduced successfully. Reduced shape: {reduced_df.shape}")

# Drop rows with missing values in critical columns
before_count = len(reduced_df)
reduced_df = reduced_df.dropna(subset=['Identifier', 'Parent', 'Classification', 'Centroid.X.um', 'Centroid.Y.um'])
after_count = len(reduced_df)
if before_count > after_count:
    logging.info(f"Dropped {before_count - after_count} rows with missing values.")

# ----------------------------------------------------------------------
# Step 3: Define cell types to analyze
# ----------------------------------------------------------------------
logging.info("Identifying cell types for interaction calculation...")

# Option 1: Extract all unique cell types
# Uncomment the next line if you want to use all unique cell types
# cell_types = sorted(reduced_df['Classification'].unique().tolist())

# Option 2: Manually specify the cell types to analyze
cell_types = [
    'M1 macrophage',
    'M1 microglia',
    'M2 macrophage',
    'M2 microglia',
    'Other cell',
    # Add more cell types as needed
]

# Log the counts of each cell type in the data for validation
cell_type_counts = reduced_df['Classification'].value_counts()
logging.info(f"Cell type distribution in the data:\n{cell_type_counts}")

logging.info(f"Selected {len(cell_types)} cell types to analyze")
logging.info(f"Will calculate {len(cell_types) * len(cell_types)} pairwise interactions")

# ----------------------------------------------------------------------
# Step 4: Calculate interaction indices grouped by BOTH Image ID and Hypoxia Class
# ----------------------------------------------------------------------
logging.info("Defining maximum interaction distance for interaction calculation...")
max_distance = 20  # Adjust as needed

# Create a combined identifier for each image-hypoxia class combination
reduced_df['Image_HypoxiaClass'] = reduced_df['Identifier'] + '_' + reduced_df['Parent'].astype(str)
logging.info(f"Created combined Image_HypoxiaClass identifier")

# Get unique combinations of Identifier and Parent
unique_image_hypoxia_combinations = reduced_df[['Identifier', 'Parent', 'Image_HypoxiaClass']].drop_duplicates()
logging.info(f"Found {len(unique_image_hypoxia_combinations)} unique image-hypoxia class combinations")

# Initialize DataFrame for interaction indices with these combined identifiers
interaction_indices_df = unique_image_hypoxia_combinations.copy()

def calculate_interaction_index_categorical(df, reference_type, target_type, max_distance=10):
    """
    Calculate interaction index between reference cells and target cells
    using a categorical cell type column.
    
    Parameters:
    -----------
    df : DataFrame
        DataFrame containing cell data with coordinates and cell types
    reference_type : str
        Value in the Classification column that identifies reference cells
    target_type : str
        Value in the Classification column that identifies target cells
    max_distance : float
        Maximum distance to consider for interaction
        
    Returns:
    --------
    float
        Interaction index value
    """
    # Select reference cells and target cells based on the Classification column
    reference_cells = df[df['Classification'] == reference_type]
    target_cells = df[df['Classification'] == target_type]
    
    if len(reference_cells) == 0 or len(target_cells) == 0:
        return np.nan
    
    try:
        reference_coords = reference_cells[['Centroid.X.um', 'Centroid.Y.um']].values
        target_coords = target_cells[['Centroid.X.um', 'Centroid.Y.um']].values
        
        reference_tree = cKDTree(reference_coords)
        target_tree = cKDTree(target_coords)
        
        # Find, for each reference cell, indices of target cells within max_distance
        interaction_counts = reference_tree.query_ball_tree(target_tree, r=max_distance)
        
        # Count reference cells with at least one nearby target cell
        reference_cells_with_interaction = sum(len(counts) > 0 for counts in interaction_counts)
        
        total_reference_cells = len(reference_cells)
        total_target_cells = len(target_cells)
        
        if total_reference_cells * total_target_cells > 0:
            interaction_index = reference_cells_with_interaction / np.sqrt(total_reference_cells * total_target_cells)
            return interaction_index
        else:
            return 0
    except Exception as e:
        logging.error(f"Error in interaction calculation: {str(e)}")
        return np.nan

# Calculate interaction indices for ALL pairwise combinations
total_combinations = len(cell_types) * len(cell_types)
current_combination = 0

for reference_type in cell_types:
    for target_type in cell_types:
        current_combination += 1
        logging.info(f"Calculating interaction index for {reference_type} -> {target_type} ({current_combination}/{total_combinations})...")
        
        # Create an empty list to store results
        interaction_results = []
        
        # Get unique combinations of Identifier and Parent
        unique_combinations = reduced_df[['Identifier', 'Parent']].drop_duplicates()
        
        # Loop through each unique combination and calculate the interaction index
        for _, row in unique_combinations.iterrows():
            identifier = row['Identifier']
            parent = row['Parent']
            
            # Filter the dataframe for this specific image and hypoxia class
            subset_df = reduced_df[(reduced_df['Identifier'] == identifier) & 
                                  (reduced_df['Parent'] == parent)]
            
            # Calculate the interaction index
            index_value = calculate_interaction_index_categorical(
                subset_df, 
                reference_type, 
                target_type, 
                max_distance
            )
            
            # Store the result
            interaction_results.append({
                'Identifier': identifier,
                'Parent': parent,
                f'{target_type}_{reference_type}_Interaction_Index': index_value
            })
        
        # Convert results to a DataFrame
        interaction_index_df = pd.DataFrame(interaction_results)
        
        # Merge with the main interaction indices DataFrame
        interaction_indices_df = pd.merge(
            interaction_indices_df,
            interaction_index_df,
            on=['Identifier', 'Parent'],
            how='left'
        )
        
        logging.info(f"Completed calculation for {reference_type} -> {target_type}.")

# Check the results for missing values
missing_values = interaction_indices_df.isna().sum()
logging.info(f"Missing values in interaction indices DataFrame:\n{missing_values}")

# ----------------------------------------------------------------------
# Step 5: Save the final DataFrame with interaction indices
# ----------------------------------------------------------------------
final_df = interaction_indices_df

output_csv_path = "/scratch/svc_td_cri/projects/multiplex/Old_narvi/hypoxia/intensities/interactions/different_distances/cell_interactions_all_pairwise_by_hypoxia_20.csv"
final_df.to_csv(output_csv_path, index=False)
logging.info(f"Final DataFrame saved to CSV at: {output_csv_path}")

# ----------------------------------------------------------------------
# Step 6: Generate summary statistics by hypoxia class
# ----------------------------------------------------------------------
logging.info("Generating summary statistics by hypoxia class...")

# Create a summary DataFrame with mean interaction indices by hypoxia class
summary_columns = [col for col in final_df.columns if 'Interaction_Index' in col]
if summary_columns:
    hypoxia_summary = final_df.groupby('Parent')[summary_columns].mean().reset_index()
    
    summary_csv_path = "/scratch/svc_td_cri/projects/multiplex/Old_narvi/hypoxia/intensities/interactions/different_distances/interaction_summary_all_pairwise_by_hypoxia_20.csv"
    hypoxia_summary.to_csv(summary_csv_path, index=False)
    logging.info(f"Hypoxia class summary saved to CSV at: {summary_csv_path}")
    
    # Display summary statistics in the log
    logging.info("Summary of interactions by hypoxia class (showing first 10 interactions):")
    for i, col in enumerate(summary_columns[:10]):
        try:
            class_means = final_df.groupby('Parent')[col].mean()
            logging.info(f"\n{col} average by hypoxia class:\n{class_means}")
        except Exception as e:
            logging.error(f"Error generating summary for {col}: {str(e)}")
    
    if len(summary_columns) > 10:
        logging.info(f"\n... and {len(summary_columns) - 10} more interaction indices")

# Show a summary of the final results
logging.info(f"Calculation complete. Generated interaction indices for {len(summary_columns)} pairwise combinations.")
logging.info(f"Final DataFrame shape: {final_df.shape}")
logging.info(f"Number of interaction index columns: {len(summary_columns)}")
logging.info("Script execution completed successfully.")
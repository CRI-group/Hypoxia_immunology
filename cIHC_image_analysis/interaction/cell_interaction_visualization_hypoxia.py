import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages
from scipy.spatial import cKDTree
import matplotlib.colors as mcolors
import matplotlib.patches as mpatches
import logging

# Set up logging
logging.basicConfig(level=logging.INFO, format='%(asctime)s - %(levelname)s - %(message)s')

# -----------------------------------------------------------------------------
# Configuration: Define cell types and their visualization properties
# -----------------------------------------------------------------------------
CELL_TYPE_CONFIG = {
    'M1 macrophage': {
        'color': 'red',
        'label': 'M1 macrophage',
        'marker_normal': 'o',
        'marker_interacting': '*',
        'size_normal': 30,
        'size_interacting': 80
    },
    'M2 macrophage': {
        'color': 'green',
        'label': 'M2 macrophage', 
        'marker_normal': 'o',
        'marker_interacting': '*',
        'size_normal': 30,
        'size_interacting': 80
    },
    'M1 microglia': {
        'color': 'orange',
        'label': 'M1 microglia',
        'marker_normal': 'o', 
        'marker_interacting': '*',
        'size_normal': 30,
        'size_interacting': 80
    },
    'M2 microglia': {
        'color': 'purple',
        'label': 'M2 microglia',
        'marker_normal': 'o',
        'marker_interacting': '*', 
        'size_normal': 30,
        'size_interacting': 80
    },
    'Other cell': {
        'color': 'blue', 
        'label': 'Other cell',
        'marker_normal': 'o',
        'marker_interacting': '*',
        'size_normal': 30,
        'size_interacting': 80
    },
    # Add more cell types here as needed with their visualization properties
    # Template for new cell types:
    # 'New Cell Type': {
    #     'color': 'color_name',
    #     'label': 'Display Name',
    #     'marker_normal': 'o',
    #     'marker_interacting': '*',
    #     'size_normal': 30,
    #     'size_interacting': 80
    # }
}

# Define which cell types to visualize and analyze
# You can easily modify this list to change what gets plotted
CELL_TYPES_TO_PLOT = ['M1 macrophage', 'Other cell']  # Change these to any cell types you want

# Define the reference cell type for interaction calculation
REFERENCE_CELL_TYPE = 'Other cell'  # Change this to any cell type you want as reference

# -----------------------------------------------------------------------------
# Step 1: Load the original cell data and interaction data
# -----------------------------------------------------------------------------
def load_data():
    # Load the original data with all cells
    input_df_path = "/lustre/compbio/projects/Multiplex_ihc_image_analysis/hypoxia/intensities/combined_Hypoxia_alldata.csv"
    df = pd.read_csv(input_df_path, low_memory=False)
    logging.info(f"Original data loaded. Shape: {df.shape}")
    
    # Load the interaction data
    interaction_df_path = "/lustre/compbio/projects/Multiplex_ihc_image_analysis/hypoxia/intensities/interactions/cell_interactions_by_hypoxia_class_OtherCells_20.csv"
    interactions = pd.read_csv(interaction_df_path)
    logging.info(f"Interaction data loaded. Shape: {interactions.shape}")
    
    return df, interactions

# -----------------------------------------------------------------------------
# Utility function: Get available cell types from the data
# -----------------------------------------------------------------------------
def get_available_cell_types(df):
    """
    Get all unique cell types available in the dataset
    
    Parameters:
    -----------
    df : DataFrame
        DataFrame containing all cell data
        
    Returns:
    --------
    list
        List of unique cell types
    """
    cell_types = sorted(df['Classification'].unique())
    logging.info(f"Available cell types in the dataset: {cell_types}")
    return cell_types

# -----------------------------------------------------------------------------
# Utility function: Validate cell types
# -----------------------------------------------------------------------------
def validate_cell_types(df, cell_types_to_plot, reference_cell_type):
    """
    Validate that the specified cell types exist in the data
    
    Parameters:
    -----------
    df : DataFrame
        DataFrame containing all cell data
    cell_types_to_plot : list
        List of cell types to validate
    reference_cell_type : str
        Reference cell type to validate
        
    Returns:
    --------
    bool
        True if all cell types are valid, False otherwise
    """
    available_cell_types = get_available_cell_types(df)
    
    # Check if all cell types to plot are available
    invalid_types = [ct for ct in cell_types_to_plot if ct not in available_cell_types]
    if invalid_types:
        logging.error(f"Invalid cell types specified: {invalid_types}")
        logging.error(f"Available cell types: {available_cell_types}")
        return False
    
    # Check if reference cell type is available
    if reference_cell_type not in available_cell_types:
        logging.error(f"Invalid reference cell type: {reference_cell_type}")
        logging.error(f"Available cell types: {available_cell_types}")
        return False
    
    return True

# -----------------------------------------------------------------------------
# Step 2: Identify interacting cells for a specific image with hypoxia information
# -----------------------------------------------------------------------------
def identify_interacting_cells(df, image_id, cell_types_to_analyze, reference_cell_type, max_distance=10):
    """
    For a specific image, identify which cells are interacting
    
    Parameters:
    -----------
    df : DataFrame
        DataFrame containing all cell data
    image_id : str
        The identifier of the image to analyze
    cell_types_to_analyze : list
        List of cell types to analyze for interactions
    reference_cell_type : str
        The reference cell type for interaction calculation
    max_distance : float
        The maximum distance to consider for interaction
        
    Returns:
    --------
    DataFrame
        Original DataFrame with an additional column indicating interacting cells
    """
    # Filter for the specific image
    image_df = df[df['Identifier'] == image_id].copy()
    
    if len(image_df) == 0:
        logging.error(f"No data found for image ID: {image_id}")
        return None
    
    # Get unique hypoxia classes in this image
    hypoxia_classes = image_df['Parent'].unique()
    logging.info(f"Found {len(hypoxia_classes)} hypoxia classes in image {image_id}: {hypoxia_classes}")
    
    # Initialize a column for interaction status (0 = not interacting, 1 = interacting)
    image_df['Is_Interacting'] = 0
    
    # Process each hypoxia class separately
    for hypoxia_class in hypoxia_classes:
        # Filter for this hypoxia class
        hypoxia_df = image_df[image_df['Parent'] == hypoxia_class]
        
        # Get reference cells
        reference_cells = hypoxia_df[hypoxia_df['Classification'] == reference_cell_type]
        
        # Process each target cell type
        for target_cell_type in cell_types_to_analyze:
            if target_cell_type == reference_cell_type:
                continue  # Skip self-interaction
                
            target_cells = hypoxia_df[hypoxia_df['Classification'] == target_cell_type]
            
            logging.info(f"Hypoxia class {hypoxia_class}: Found {len(reference_cells)} {reference_cell_type} and {len(target_cells)} {target_cell_type}")
            
            if len(reference_cells) > 0 and len(target_cells) > 0:
                # Create KD-Trees for spatial queries
                ref_coords = reference_cells[['Centroid.X.um', 'Centroid.Y.um']].values
                target_coords = target_cells[['Centroid.X.um', 'Centroid.Y.um']].values
                
                ref_tree = cKDTree(ref_coords)
                target_tree = cKDTree(target_coords)
                
                # Find interactions in both directions
                # Reference cells with nearby target cells
                ref_interactions = ref_tree.query_ball_tree(target_tree, r=max_distance)
                
                # Target cells with nearby reference cells
                target_interactions = target_tree.query_ball_tree(ref_tree, r=max_distance)
                
                # Mark reference cells that interact with target cells
                for i, interactions in enumerate(ref_interactions):
                    if len(interactions) > 0:
                        cell_index = reference_cells.iloc[i].name
                        image_df.at[cell_index, 'Is_Interacting'] = 1
                
                # Mark target cells that interact with reference cells
                for i, interactions in enumerate(target_interactions):
                    if len(interactions) > 0:
                        cell_index = target_cells.iloc[i].name
                        image_df.at[cell_index, 'Is_Interacting'] = 1
                        
                logging.info(f"Identified interactions between {reference_cell_type} and {target_cell_type} in hypoxia class {hypoxia_class}")
    
    total_interacting = image_df['Is_Interacting'].sum()
    logging.info(f"Total interacting cells identified: {total_interacting}")
    
    return image_df

# -----------------------------------------------------------------------------
# Step 3: Visualize cells in the image with hypoxia information
# -----------------------------------------------------------------------------
def visualize_cell_interactions(df, image_id, cell_types_to_plot, reference_cell_type, save_path=None, visualization_mode='combined'):
    """
    Create a visualization of cells in an image, highlighting interacting cells and hypoxia classes
    
    Parameters:
    -----------
    df : DataFrame
        DataFrame containing all cell data with interaction status
    image_id : str
        The identifier of the image to visualize
    cell_types_to_plot : list
        List of cell types to include in the visualization
    reference_cell_type : str
        Reference cell type for interaction calculation
    save_path : str
        Path to save the visualization
    visualization_mode : str
        'combined' - Show all hypoxia classes in one plot with different background colors
        'separate' - Create separate plots for each hypoxia class
    """
    # Validate cell types
    if not validate_cell_types(df, cell_types_to_plot, reference_cell_type):
        return None
    
    # Filter for the image and identify interacting cells
    image_df = identify_interacting_cells(df, image_id, cell_types_to_plot, reference_cell_type)
    
    if image_df is None or len(image_df) == 0:
        logging.error(f"Cannot create visualization for image {image_id}: No data")
        return
    
    # Filter for specified cell types
    target_cells = image_df[image_df['Classification'].isin(cell_types_to_plot)]
    
    # Get unique hypoxia classes
    hypoxia_classes = sorted(target_cells['Parent'].unique())
    
    # Define colors for hypoxia classes
    hypoxia_cmap = plt.cm.Blues
    hypoxia_colors = {cls: hypoxia_cmap(i / len(hypoxia_classes)) 
                      for i, cls in enumerate(hypoxia_classes)}
    
    def plot_cell_type(cell_data, cell_type, hypoxia_class=None, show_in_legend=True):
        """Helper function to plot a specific cell type"""
        if cell_type not in CELL_TYPE_CONFIG:
            logging.warning(f"Cell type {cell_type} not found in configuration. Using default settings.")
            config = {
                'color': 'gray',
                'label': cell_type,
                'marker_normal': 'o',
                'marker_interacting': '*',
                'size_normal': 30,
                'size_interacting': 80
            }
        else:
            config = CELL_TYPE_CONFIG[cell_type]
        
        # Non-interacting cells
        non_interacting = cell_data[(cell_data['Classification'] == cell_type) & 
                                   (cell_data['Is_Interacting'] == 0)]
        if len(non_interacting) > 0:
            label = f"{config['label']} (non-interacting)" if show_in_legend else ""
            plt.scatter(non_interacting['Centroid.X.um'], 
                       non_interacting['Centroid.Y.um'],
                       color=config['color'], alpha=0.2, s=config['size_normal'], 
                       marker=config['marker_normal'], label=label)
        
        # Interacting cells
        interacting = cell_data[(cell_data['Classification'] == cell_type) & 
                               (cell_data['Is_Interacting'] == 1)]
        if len(interacting) > 0:
            label = f"{config['label']} (interacting)" if show_in_legend else ""
            plt.scatter(interacting['Centroid.X.um'], 
                       interacting['Centroid.Y.um'],
                       color=config['color'], alpha=0.8, s=config['size_interacting'], 
                       marker=config['marker_interacting'], label=label)
        
        return len(non_interacting), len(interacting)
    
    if visualization_mode == 'separate':
        # Create separate plots for each hypoxia class and save to multi-page PDF
        if save_path:
            with PdfPages(save_path) as pdf:
                for hypoxia_class in hypoxia_classes:
                    # Filter for this hypoxia class
                    hypoxia_df = target_cells[target_cells['Parent'] == hypoxia_class]
                    
                    plt.figure(figsize=(10, 8))
                    
                    # Plot each cell type
                    stats_text = f"Hypoxia Class: {hypoxia_class}\n"
                    
                    for i, cell_type in enumerate(cell_types_to_plot):
                        non_int, interacting = plot_cell_type(hypoxia_df, cell_type, hypoxia_class, show_in_legend=True)
                        total = non_int + interacting
                        if total > 0:
                            stats_text += f"{cell_type}: {interacting}/{total} interacting ({interacting/total*100:.1f}%)\n"
                    
                    # Add the statistics to the plot
                    plt.text(0.02, 0.98, stats_text, transform=plt.gca().transAxes, 
                             verticalalignment='top', horizontalalignment='left',
                             bbox=dict(boxstyle='round', facecolor='white', alpha=0.8))
                    
                    # Set title and labels
                    plt.title(f'Cell Interactions in Image: {image_id} - Hypoxia Class: {hypoxia_class}', fontsize=14)
                    plt.xlabel('X Position (μm)', fontsize=12)
                    plt.ylabel('Y Position (μm)', fontsize=12)
                    plt.legend(loc='upper right')
                    plt.grid(True, alpha=0.3)
                    
                    # Add scale bar
                    x_min, x_max = plt.gca().get_xlim()
                    y_min, y_max = plt.gca().get_ylim()
                    scale_bar_length = 100  # 100μm
                    scale_bar_x = x_min + (x_max - x_min) * 0.05
                    scale_bar_y = y_min + (y_max - y_min) * 0.05
                    plt.plot([scale_bar_x, scale_bar_x + scale_bar_length], 
                             [scale_bar_y, scale_bar_y], 'k-', lw=2)
                    plt.text(scale_bar_x + scale_bar_length/2, scale_bar_y - (y_max - y_min) * 0.01, 
                             '100 μm', ha='center')
                    
                    plt.tight_layout()
                    
                    # Save the current figure to the PDF
                    pdf.savefig(plt.gcf(), dpi=300, bbox_inches='tight')
                    logging.info(f"Added hypoxia class {hypoxia_class} to PDF")
                    
                    # Close the current figure to free memory
                    plt.close()
                
                logging.info(f"All separate hypoxia visualizations saved to {save_path}")
        else:
            # If no save path, display plots normally
            for hypoxia_class in hypoxia_classes:
                # Filter for this hypoxia class
                hypoxia_df = target_cells[target_cells['Parent'] == hypoxia_class]
                
                plt.figure(figsize=(10, 8))
                
                # Plot each cell type
                stats_text = f"Hypoxia Class: {hypoxia_class}\n"
                
                for i, cell_type in enumerate(cell_types_to_plot):
                    non_int, interacting = plot_cell_type(hypoxia_df, cell_type, hypoxia_class, show_in_legend=True)
                    total = non_int + interacting
                    if total > 0:
                        stats_text += f"{cell_type}: {interacting}/{total} interacting ({interacting/total*100:.1f}%)\n"
                
                # Add the statistics to the plot
                plt.text(0.02, 0.98, stats_text, transform=plt.gca().transAxes, 
                         verticalalignment='top', horizontalalignment='left',
                         bbox=dict(boxstyle='round', facecolor='white', alpha=0.8))
                
                # Set title and labels
                plt.title(f'Cell Interactions in Image: {image_id} - Hypoxia Class: {hypoxia_class}', fontsize=14)
                plt.xlabel('X Position (μm)', fontsize=12)
                plt.ylabel('Y Position (μm)', fontsize=12)
                plt.legend(loc='upper right')
                plt.grid(True, alpha=0.3)
                
                # Add scale bar
                x_min, x_max = plt.gca().get_xlim()
                y_min, y_max = plt.gca().get_ylim()
                scale_bar_length = 100  # 100μm
                scale_bar_x = x_min + (x_max - x_min) * 0.05
                scale_bar_y = y_min + (y_max - y_min) * 0.05
                plt.plot([scale_bar_x, scale_bar_x + scale_bar_length], 
                         [scale_bar_y, scale_bar_y], 'k-', lw=2)
                plt.text(scale_bar_x + scale_bar_length/2, scale_bar_y - (y_max - y_min) * 0.01, 
                         '100 μm', ha='center')
                
                plt.tight_layout()
    
    else:  # Combined visualization
        # Create a combined visualization with all hypoxia classes
        plt.figure(figsize=(14, 12))
        
        # Plot cells by hypoxia class
        stats_text = "Statistics by Hypoxia Class:\n"
        
        for i, hypoxia_class in enumerate(hypoxia_classes):
            # Get data for this hypoxia class
            hypoxia_df = target_cells[target_cells['Parent'] == hypoxia_class]
            
            # Create a background patch to represent the hypoxia area
            if len(hypoxia_df) > 0:
                min_x = hypoxia_df['Centroid.X.um'].min()
                max_x = hypoxia_df['Centroid.X.um'].max()
                min_y = hypoxia_df['Centroid.Y.um'].min()
                max_y = hypoxia_df['Centroid.Y.um'].max()
                
                # Add a slightly transparent rectangular patch for this hypoxia region
                rect = plt.Rectangle((min_x, min_y), max_x - min_x, max_y - min_y,
                                     color=hypoxia_colors[hypoxia_class], alpha=0.1, 
                                     label=f'Hypoxia Class {hypoxia_class}')
                plt.gca().add_patch(rect)
            
            # Plot each cell type for this hypoxia class
            stats_text += f"\nClass {hypoxia_class}:\n"
            
            for cell_type in cell_types_to_plot:
                non_int, interacting = plot_cell_type(hypoxia_df, cell_type, hypoxia_class, 
                                                     show_in_legend=(i == 0))  # Only show in legend for first hypoxia class
                total = non_int + interacting
                if total > 0:
                    stats_text += f"  {cell_type}: {interacting}/{total} interacting ({interacting/total*100:.1f}%)\n"
        
        # Create a custom legend for hypoxia classes
        hypoxia_patches = [mpatches.Patch(color=hypoxia_colors[cls], alpha=0.3, 
                                         label=f'Hypoxia Class {cls}') 
                          for cls in hypoxia_classes]
        
        # Add the overall statistics
        plt.text(0.02, 0.98, stats_text, transform=plt.gca().transAxes, 
                 verticalalignment='top', horizontalalignment='left',
                 bbox=dict(boxstyle='round', facecolor='white', alpha=0.8))
        
        # Set title and labels
        plt.title(f'Cell Interactions by Hypoxia Class in Image: {image_id}', fontsize=14)
        plt.xlabel('X Position (μm)', fontsize=12)
        plt.ylabel('Y Position (μm)', fontsize=12)
        
        # Combine the existing legend with hypoxia patches
        handles, labels = plt.gca().get_legend_handles_labels()
        plt.legend(handles=handles + hypoxia_patches, loc='upper right')
        
        plt.grid(True, alpha=0.3)
        
        # Add scale bar
        x_min, x_max = plt.gca().get_xlim()
        y_min, y_max = plt.gca().get_ylim()
        scale_bar_length = 100  # 100μm
        scale_bar_x = x_min + (x_max - x_min) * 0.05
        scale_bar_y = y_min + (y_max - y_min) * 0.05
        plt.plot([scale_bar_x, scale_bar_x + scale_bar_length], 
                [scale_bar_y, scale_bar_y], 'k-', lw=2)
        plt.text(scale_bar_x + scale_bar_length/2, scale_bar_y - (y_max - y_min) * 0.01, 
                '100 μm', ha='center')
        
        plt.tight_layout()
        
        # Save the figure if a save path is provided
        if save_path:
            plt.savefig(save_path, dpi=300, bbox_inches='tight')
            logging.info(f"Combined visualization saved to {save_path}")
    
    return plt.gcf()  # Return the figure

# -----------------------------------------------------------------------------
# Step 4: Create a main function to run the visualization
# -----------------------------------------------------------------------------
def main():
    # Load the data
    df, interactions = load_data()
    
    # Show available cell types
    available_cell_types = get_available_cell_types(df)
    
    # Specify the image ID to visualize
    # You can change this to whatever image ID you want to visualize
    image_id = "Astro 24_A3.tif"
    
    # *** EASILY CONFIGURABLE SECTION ***
    # Change these variables to modify what gets plotted:
    
    # 1. Define which cell types to plot (choose from available_cell_types)
    cell_types_to_plot = ['M1 macrophage', 'Other cell']  # Fixed: was 'Other cells'
    
    # 2. Define the reference cell type for interaction calculation
    reference_cell_type = 'Other cell'  # Fixed: was 'Other cells'
    
    # Alternative: You can also specify them directly here:
    # cell_types_to_plot = ['M1 macrophage', 'M2 macrophage', 'M1 microglia']
    # reference_cell_type = 'M1 macrophage'
    
    logging.info(f"Plotting cell types: {cell_types_to_plot}")
    logging.info(f"Using reference cell type: {reference_cell_type}")
    
    # Create both types of visualizations
    
    # 1. Combined visualization (all hypoxia classes in one plot)
    combined_save_path = "/lustre/compbio/projects/Multiplex_ihc_image_analysis/hypoxia/intensities/figures/interaction/cell_interaction_Hypoxia_combined_OtherCell_M1macros_20.pdf"
    
    logging.info("Creating combined visualization...")
    visualize_cell_interactions(df, image_id, cell_types_to_plot, reference_cell_type, 
                               combined_save_path, visualization_mode='combined')
    
    # 2. Separate visualizations (one plot per hypoxia class in a single multi-page PDF)
    separate_save_path = "/lustre/compbio/projects/Multiplex_ihc_image_analysis/hypoxia/intensities/figures/interaction/cell_interaction_Hypoxia_separate_OtherCell_M1macros_20.pdf"
    
    logging.info("Creating separate visualizations for each hypoxia class in a single PDF...")
    visualize_cell_interactions(df, image_id, cell_types_to_plot, reference_cell_type, 
                               separate_save_path, visualization_mode='separate')
    
    logging.info("Visualization process completed.")

# -----------------------------------------------------------------------------
# Step 5: Additional utility function to create custom visualizations
# -----------------------------------------------------------------------------
def create_custom_visualization(df, image_id, cell_types_to_plot, reference_cell_type, 
                               save_path=None, visualization_mode='combined'):
    """
    Create a custom visualization with specified parameters
    
    Parameters:
    -----------
    df : DataFrame
        DataFrame containing all cell data
    image_id : str
        The identifier of the image to visualize
    cell_types_to_plot : list
        List of cell types to include in the visualization
    reference_cell_type : str
        Reference cell type for interaction calculation
    save_path : str
        Path to save the visualization
    visualization_mode : str
        'combined' or 'separate'
    
    Example usage:
    --------------
    # Plot M1 and M2 macrophages with M1 as reference
    create_custom_visualization(df, "Astro 24_A3.tif", 
                               ['M1 macrophage', 'M2 macrophage'], 
                               'M1 macrophage', 
                               'custom_viz.png', 'combined')
    """
    logging.info(f"Creating custom visualization for cell types: {cell_types_to_plot}")
    logging.info(f"Reference cell type: {reference_cell_type}")
    
    return visualize_cell_interactions(df, image_id, cell_types_to_plot, reference_cell_type, 
                                     save_path, visualization_mode)

# -----------------------------------------------------------------------------
# Step 6: Additional utility function to explore hypoxia classes
# -----------------------------------------------------------------------------
def explore_hypoxia_classes(df, image_id=None):
    """
    Explore the hypoxia classes available in the data
    
    Parameters:
    -----------
    df : DataFrame
        DataFrame containing all cell data
    image_id : str, optional
        Specific image to explore, if None explores all data
    """
    if image_id:
        data_subset = df[df['Identifier'] == image_id]
        logging.info(f"Exploring hypoxia classes for image: {image_id}")
    else:
        data_subset = df
        logging.info("Exploring hypoxia classes for all data")
    
    # Get unique hypoxia classes
    hypoxia_classes = data_subset['Parent'].unique()
    logging.info(f"Found {len(hypoxia_classes)} unique hypoxia classes: {sorted(hypoxia_classes)}")
    
    # Count cells by hypoxia class and cell type
    hypoxia_cell_counts = data_subset.groupby(['Parent', 'Classification']).size().unstack(fill_value=0)
    logging.info(f"Cell counts by hypoxia class and cell type:\n{hypoxia_cell_counts}")
    
    return hypoxia_classes, hypoxia_cell_counts

# Run the main function if this script is executed
if __name__ == "__main__":
    # First explore the data to understand hypoxia classes
    df, interactions = load_data()
    explore_hypoxia_classes(df, "Astro 24_A3.tif")
    
    # Then run the main visualization
    main()
    
    # Example of how to create custom visualizations:
    logging.info("\n" + "="*50)
    logging.info("EXAMPLE: Creating custom visualizations")
    logging.info("="*50)
    
    # Example 1: Plot M1 and M2 macrophages
    # create_custom_visualization(df, "Allglioma_B8.tif", 
    #                            ['M1 macrophage', 'M2 macrophage'], 
    #                            'M1 macrophage', 
    #                            '/path/to/m1_m2_viz.png', 'combined')
    
    # Example 2: Plot all microglia types
    # create_custom_visualization(df, "Allglioma_B8.tif", 
    #                            ['M1 microglia', 'M2 microglia'], 
    #                            'M1 microglia', 
    #                            '/path/to/microglia_viz.png', 'separate')
    
    # Example 3: Plot any combination of cell types
    # create_custom_visualization(df, "Allglioma_B8.tif", 
    #                            ['M1 macrophage', 'M1 microglia', 'M2 macrophage'], 
    #                            'M2 macrophage', 
    #                            '/path/to/custom_viz.png', 'combined')
    
    logging.info("Uncomment the example lines above to create custom visualizations!")

# =============================================================================
# QUICK CONFIGURATION GUIDE:
# =============================================================================
# 
# To easily change what gets plotted, modify these variables at the top:
# 
# 1. CELL_TYPES_TO_PLOT: List of cell types to visualize
#
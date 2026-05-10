#!/usr/bin/env python3
# Generates clustering visualizations from APD distance matrices
# Parses distance CSVs and constructs hierarchical dendrograms and heatmaps

import csv
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import seaborn as sns
from scipy.cluster import hierarchy
import scipy.spatial.distance as ssd

def plot_matrices(csv_path, summary_figs_dir):
    """Processes a distance matrix CSV and generates a hierarchical dendrogram and heatmap"""
    print(f"\n  -> Processing {csv_path.name}...")
    
    # Reads the distance matrix
    with open(csv_path, "r") as f:
        reader = csv.reader(f)
        header = next(reader)
        structure_names = header[1:]
        
        matrix = []
        for row in reader:
            matrix.append([float(x) for x in row[1:]])
            
    matrix = np.array(matrix)
    

    # Converts redundant square matrix to a condensed form required by SciPy
    condensed_dist = ssd.squareform(matrix)
    
    # Performs UPGMA clustering (average linkage)
    Z = hierarchy.linkage(condensed_dist, method='average')
    
    fig, ax = plt.subplots(figsize=(12, 8))
    
    # Renders the dendrogram
    # orientation='left' puts labels on the right side for easier reading
    hierarchy.dendrogram(
        Z, 
        labels=structure_names, 
        orientation='left', 
        leaf_font_size=12,
        link_color_func=lambda k: 'black',
        ax=ax
    )
    
    # Removes the top, right, and left spines for a cleaner aesthetic
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)
    ax.spines['left'].set_visible(False)
    
    plt.title(f"Hierarchical Clustering Dendrogram ({csv_path.stem})", fontsize=14)
    plt.xlabel("Amyloid Packing Difference (%)", fontsize=12)
    plt.tight_layout()
    
    dendro_path = summary_figs_dir / f"{csv_path.stem}_dendrogram.svg"
    fig.savefig(dendro_path, dpi=300)
    plt.close(fig)
    print(f"    -> Generated Dendrogram: {dendro_path.name}")

    # Reorders the matrix and labels based on the hierarchical clustering
    leaves = hierarchy.leaves_list(Z)
    reordered_matrix = matrix[leaves, :][:, leaves]
    reordered_names = [structure_names[i] for i in leaves]

    # Generates a distance matrix heatmap
    plt.figure(figsize=(10, 8))
    # Uses a yellow-orange-red colormap which fits distance metrics well
    sns.heatmap(reordered_matrix, xticklabels=reordered_names, yticklabels=reordered_names, 
                cmap="YlOrRd", cbar_kws={'label': 'APD (%)'})
    plt.title(f"APD Distance Matrix Heatmap ({csv_path.stem})", fontsize=14)
    plt.tight_layout()
    
    heatmap_path = summary_figs_dir / f"{csv_path.stem}_heatmap.svg"
    plt.savefig(heatmap_path, dpi=300)
    plt.close()
    print(f"    -> Generated Heatmap: {heatmap_path.name}")

def main():
    """Executes the clustering and plotting workflow"""
    comparisons_dir = Path("Comparisons")
    summary_figs_dir = comparisons_dir / "Summary_Figures"
    
    if not comparisons_dir.exists():
        print(f"  -> [!] Error: Directory '{comparisons_dir}' not found")
        return
        
    summary_figs_dir.mkdir(parents=True, exist_ok=True)
    
    # Processes both the XY and XYZ matrices automatically
    csv_files = list(comparisons_dir.glob("apd_distance_matrix_*.csv"))
    if not csv_files:
        print(f"  -> [!] No distance matrix CSVs found in '{comparisons_dir}'")
        return
        
    for csv_file in csv_files:
        plot_matrices(csv_file, summary_figs_dir)
        
    print("\n  -> All plots generated successfully! Check the Comparisons/Summary_Figures directory\n")

if __name__ == "__main__":
    main()
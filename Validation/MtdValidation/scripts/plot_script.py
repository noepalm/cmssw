#!/usr/bin/env python3
"""
Plot script for MTD MergedCluster validation histograms
Uses mplhep and uproot to read ROOT files and create publication-quality plots
"""

import uproot
import matplotlib.pyplot as plt
import mplhep as hep
import numpy as np
import os
import argparse
from pathlib import Path

# Set matplotlib style
plt.style.use(hep.style.CMS)

def plot_histogram(hist, name, output_dir, log_scale=False):
    """
    Plot a single histogram with CMS style
    """
    # strip name of complementary bits
    stripped_name = name.replace("h_", "")

    fig, ax = plt.subplots(figsize=(10, 8))
    
    # Get bin centers and values
    bin_edges = hist.axis().edges()
    bin_centers = hist.axis().centers()
    values = hist.values()
    errors = np.sqrt(hist.variances()) if hasattr(hist, 'variances') else None
    
    # Create the plot
    hep.histplot(values, bins=bin_edges, yerr=errors, ax=ax, 
                histtype='step', linewidth=2, label='Data')
    
    # Set labels and title
    ax.set_xlabel(stripped_name)
    # ax.set_xlabel(hist.axis().label if hasattr(hist.axis(), 'label') else 'Value')
    ax.set_ylabel('Events')
    ax.set_title(stripped_name.replace('_', ' ').title(), pad=40)
    
    # Set log scale if requested
    if log_scale:
        ax.set_yscale('log')
    
    # Add CMS label
    hep.cms.label("Simulation", data=False, ax=ax)
    
    # Add statistics box
    entries = int(np.sum(values))
    mean = np.average(bin_centers, weights=values) if entries > 0 else 0
    ax.text(0.75, 0.85, f'Entries: {entries}\nMean: {mean:.3f}', 
            transform=ax.transAxes, fontsize=10)
            # bbox=dict(boxstyle="round,pad=0.3", facecolor="white", alpha=0.8))
    
    # Save the plot
    log_suffix = "_log" if log_scale else ""
    output_path = output_dir / f"{name}{log_suffix}.png"
    plt.tight_layout()
    plt.savefig(output_path, dpi=300, bbox_inches='tight')
    plt.close()
    
    print(f"Saved plot: {output_path}")

def plot_2d_histogram(hist, name, output_dir):
    """
    Plot a 2D histogram with CMS style
    """
    # strip name of complementary bits
    stripped_name = name.replace("h_", "").replace("_2D", "")


    fig, ax = plt.subplots(figsize=(12, 10))
    
    # Get histogram data
    values = hist.values()
    x_edges = hist.axes[0].edges()
    y_edges = hist.axes[1].edges()
    
    # Create 2D plot
    im = ax.imshow(values.T, origin='lower', aspect='auto',
                   extent=[x_edges[0], x_edges[-1], y_edges[0], y_edges[-1]],
                   cmap='viridis')
    
    # Add colorbar
    cbar = plt.colorbar(im, ax=ax)
    cbar.set_label('Events')
    
    # Set labels and title
    x_label = stripped_name.split('_vs_')[1].replace('_', ' ') if '_vs_' in stripped_name else 'x'
    y_label = stripped_name.split('_vs_')[0].replace('_', ' ') if '_vs_' in stripped_name else 'y'
    ax.set_xlabel(x_label)
    ax.set_ylabel(y_label)
    # ax.set_xlabel(hist.axes[0].label if hasattr(hist.axes[0], 'label') else 'X')
    # ax.set_ylabel(hist.axes[1].label if hasattr(hist.axes[1], 'label') else 'Y')
    ax.set_title(stripped_name.replace('_', ' ').title(), pad=40)
    
    # Add CMS label
    hep.cms.label("Simulation", data=False, ax=ax)
    
    # Save the plot
    output_path = output_dir / f"{name}_2D.png"
    plt.tight_layout()
    plt.savefig(output_path, dpi=300, bbox_inches='tight')
    plt.close()
    
    print(f"Saved 2D plot: {output_path}")

def main():
    parser = argparse.ArgumentParser(description='Plot MTD MergedCluster validation histograms')
    parser.add_argument('input_file', nargs='?', default='mtd_mergedcluster_validation.root',
                       help='Input ROOT file (default: mtd_mergedcluster_validation.root)')
    parser.add_argument('--output-dir', '-o', default='/eos/home-n/npalmeri/www/MTD/MergedCluster/PR_plots/vali_hists',
                       help='Output directory for plots (default: plots)')
    parser.add_argument('--log-scale', '-l', action='store_true',
                       help='Use log scale for y-axis')
    parser.add_argument('--directory', '-d', default='mtdMergedClusterValidation',
                       help='Directory in ROOT file to read histograms from')
    
    args = parser.parse_args()
    
    # Check if input file exists
    if not os.path.exists(args.input_file):
        print(f"Error: Input file {args.input_file} not found!")
        return 1
    
    # Create output directory
    output_dir = Path(args.output_dir)
    output_dir.mkdir(exist_ok=True)
    
    print(f"Reading histograms from {args.input_file}...")
    print(f"Looking in directory: {args.directory}")
    
    try:
        # Open ROOT file
        with uproot.open(args.input_file) as file:
            # Check if the directory exists
            if args.directory not in file:
                print(f"Error: Directory '{args.directory}' not found in ROOT file!")
                print("Available directories:")
                for key in file.keys():
                    print(f"  {key}")
                return 1
            
            # Navigate to the validation directory
            validation_dir = file[args.directory]
            
            # Get all histogram keys
            hist_keys = [key for key in validation_dir.keys() if validation_dir[key].classname.startswith('TH')]
            
            if not hist_keys:
                print(f"No histograms found in directory '{args.directory}'!")
                return 1
            
            print(f"Found {len(hist_keys)} histograms:")
            for key in hist_keys:
                print(f"  {key}")
            
            # Plot each histogram
            for hist_key in hist_keys:
                try:
                    hist = validation_dir[hist_key]
                    hist_name = hist_key.split(';')[0]  # Remove cycle number
                    
                    print(f"\nProcessing histogram: {hist_name}")
                    
                    # Check if it's a 2D histogram
                    if hist.classname in ['TH2D', 'TH2F', 'TH2I']:
                        plot_2d_histogram(hist, hist_name, output_dir)
                    else:
                        plot_histogram(hist, hist_name, output_dir, args.log_scale)
                        
                except Exception as e:
                    print(f"Error plotting {hist_key}: {e}")
                    continue
    
    except Exception as e:
        print(f"Error reading ROOT file: {e}")
        return 1
    
    print(f"\nAll plots saved to: {output_dir}")
    print("Done!")
    
    return 0

if __name__ == "__main__":
    exit(main())
#!/usr/bin/env python3
"""
Plot script for MTD SuperCluster validation from TTree
Uses mplhep and uproot to read ROOT TTree and create publication-quality plots
"""

import uproot
import matplotlib.pyplot as plt
import mplhep as hep
import numpy as np
import os
import argparse
from pathlib import Path
from importlib import import_module

# Set matplotlib style
plt.style.use(hep.style.CMS)

def apply_cuts(tree, cuts, branches_needed):
    """
    Apply cuts to the tree data and return filtered arrays
    
    Args:
        tree: uproot TTree object
        cuts: dict with branch names as keys and lambda functions as values
        branches_needed: list of branch names that will be used for plotting
    
    Returns:
        dict of filtered numpy arrays for each branch
    """
    if cuts is None:
        # No cuts, return all data
        result = {}
        for branch in branches_needed:
            if branch == 'sc_n':
                result[branch] = tree[branch].array(library='np')
            else:
                data = tree[branch].array(library='ak')
                result[branch] = np.concatenate([event_data for event_data in data if len(event_data) > 0])
        return result
    
    # Read all needed branches including those used for cuts
    all_branches = set(branches_needed) | set(cuts.keys())
    branch_data = {}
    
    for branch in all_branches:
        if branch == 'sc_n':
            branch_data[branch] = tree[branch].array(library='np')
        else:
            data = tree[branch].array(library='ak')
            branch_data[branch] = np.concatenate([event_data for event_data in data if len(event_data) > 0])
    
    # Apply cuts sequentially - ensure all arrays are numpy arrays
    first_supercluster_branch = None
    for branch in all_branches:
        if branch != 'sc_n':
            first_supercluster_branch = branch
            break
    
    if first_supercluster_branch is None:
        # Only sc_n branch, return as-is
        result = {}
        for branch in branches_needed:
            result[branch] = branch_data[branch]
        return result
    
    mask = np.ones(len(branch_data[first_supercluster_branch]), dtype=bool)
    
    for cut_branch, cut_func in cuts.items():
        if cut_branch in branch_data and cut_branch != 'sc_n':
            # Convert to numpy array and apply cut
            cut_data = np.array(branch_data[cut_branch])
            cut_result = cut_func(cut_data)
            # Ensure cut_result is a numpy array
            if hasattr(cut_result, 'to_numpy'):
                cut_result = cut_result.to_numpy()
            else:
                cut_result = np.array(cut_result)
            mask = mask & cut_result
    
    # Apply mask to all branches except sc_n (which is event-level)
    result = {}
    for branch in branches_needed:
        if branch == 'sc_n':
            result[branch] = branch_data[branch]
        else:
            result[branch] = branch_data[branch][mask]
    
    return result

def plot_1d_from_tree(tree, plot_def, output_dir):
    """
    Create and plot a 1D histogram from TTree data
    """
    name = plot_def['name']
    branch = plot_def['branches'][0]
    
    fig, ax = plt.subplots(figsize=(10, 8))
    
    # Apply cuts and get filtered data
    cuts = plot_def.get('cuts')
    branches_needed = plot_def['branches'][:]
    
    # Check if type breakdown is requested and we have cluster-level data
    type_breakdown = plot_def.get('type_breakdown', False)
    if type_breakdown:
        # Add sc_clusterType to branches needed
        if 'sc_clusterType' not in branches_needed:
            branches_needed.append('sc_clusterType')
    
    data_dict = apply_cuts(tree, cuts, branches_needed)
    data = data_dict[branch]

    # Flatten data
    data = np.ravel(data)
    
    # Create histogram
    if 'binning' in plot_def and plot_def['binning'] is not None:
        bins = plot_def['binning']
    else:
        bins = 'auto'  # Let matplotlib choose
    
    # Create histogram and get bin contents for error calculation
    counts, bin_edges = np.histogram(data, bins=bins)
    bin_centers = (bin_edges[:-1] + bin_edges[1:]) / 2
    
    # Calculate under/overflow counts
    underflow = np.sum(data < bin_edges[0])
    overflow = np.sum(data > bin_edges[-1])
    
    # Calculate statistical uncertainties (sqrt of counts)
    uncertainties = np.sqrt(counts)
    
    # Determine colors based on whether type breakdown is requested
    default_color = 'black' if type_breakdown else '#5790fc'
    
    # Plot histogram with error bars (total/cumulative)
    ax.errorbar(bin_centers, counts, yerr=uncertainties, fmt='o', 
                markersize=3, capsize=2, capthick=1, 
                linewidth=1.5, alpha=0.8, label='Total', color=default_color)
    
    # Also plot the histogram outline for better visibility
    ax.hist(data, bins=bins, histtype='step', linewidth=2, 
            alpha=0.6, color=default_color, label='_nolegend_')
    
    # Add type breakdown if requested
    if type_breakdown and 'sc_clusterType' in data_dict:
        cluster_types = np.ravel(data_dict['sc_clusterType'])
        
        # Ensure same length as data
        min_len = min(len(data), len(cluster_types))
        data_for_types = data[:min_len]
        cluster_types = cluster_types[:min_len]
        
        # Define cluster type names and colors
        cluster_type_names = {0: 'Primary', 1: 'Secondary', 2: 'Looper', 3: 'Calo backscatter'}
        cluster_colors = ['#e42536', '#5790fc', '#7a21dd', '#9c9ca1']  # Red, Blue, Purple, Gray
        
        # Plot each cluster type separately
        for cluster_type in sorted(set(cluster_types)):
            if cluster_type in cluster_type_names:
                type_mask = cluster_types == cluster_type
                type_data = data_for_types[type_mask]
                
                if len(type_data) > 0:
                    color = cluster_colors[cluster_type % len(cluster_colors)]
                    type_name = cluster_type_names[cluster_type]
                    
                    # Create histogram for this type
                    type_counts, _ = np.histogram(type_data, bins=bins)
                    type_uncertainties = np.sqrt(type_counts)
                    
                    # Plot with error bars
                    ax.errorbar(bin_centers, type_counts, yerr=type_uncertainties, 
                               fmt='s', markersize=2, capsize=1, capthick=0.5,
                               linewidth=1, alpha=0.7, label=f'{type_name} ({len(type_data)})', 
                               color=color)
                    
                    # Also plot step histogram
                    ax.hist(type_data, bins=bins, histtype='step', linewidth=1.5, 
                           alpha=0.5, color=color, label='_nolegend_')
    
    # Set labels and title
    ax.set_xlabel(plot_def['xlabel'])
    ax.set_ylabel(plot_def['ylabel'])
    # ax.set_title(plot_def['xlabel'].replace('[', '(').replace(']', ')'), pad=40)

    # Set log scale if requested
    if plot_def.get('logy_scale', False):
        ax.set_yscale('log')
        # Ensure we don't have zero counts for log scale
        # For error bars, we need to handle zero counts differently
        nonzero_mask = counts > 0
        if np.any(nonzero_mask):
            ax.set_ylim(bottom=0.1)
        else:
            ax.set_ylim(bottom=0.01)  # If all bins are zero
    
    if plot_def.get('logx_scale', False):
        ax.set_xscale('log')
    
    # Add CMS label
    hep.cms.label("Simulation", data=False,  ax=ax)
    
    # Add legend if type breakdown is shown
    if type_breakdown and 'sc_clusterType' in data_dict:
        ax.legend(fontsize=9)
    
    # Add statistics
    entries = len(data)
    mean_val = np.mean(data) if len(data) > 0 else 0
    std_val = np.std(data) if len(data) > 0 else 0
    
    stats_text = f'Entries: {entries}\nMean: {mean_val:.3f}\nStd: {std_val:.3f}'
    stats_text += f'\nUnderflow: {underflow}\nOverflow: {overflow}'
    if cuts is not None:
        stats_text += f'\n(with cuts applied)'
    
    # Position stats box based on whether legend is present
    stats_x = 0.75 if not (type_breakdown and 'sc_clusterType' in data_dict) else 0.02
    stats_y = 0.85 if not (type_breakdown and 'sc_clusterType' in data_dict) else 0.98
    v_align = 'bottom' if not (type_breakdown and 'sc_clusterType' in data_dict) else 'top'
    
    ax.text(stats_x, stats_y, stats_text, transform=ax.transAxes, fontsize=10,
            verticalalignment=v_align,
            bbox=dict(boxstyle="round,pad=0.3", facecolor="white", alpha=0.8))
    
    # Save the plot
    # cuts_suffix = "_cut" if cuts is not None else ""
    output_path = output_dir / f"{name}.png"
    plt.tight_layout()
    plt.savefig(output_path, dpi=300, bbox_inches='tight')
    plt.close()
    
    print(f"Saved 1D plot: {output_path}")


def plot_2d_from_tree(tree, plot_def, output_dir, logy_scale=False):
    """
    Create and plot a 2D histogram from TTree data
    """
    name = plot_def['name']
    x_branch, y_branch = plot_def['branches']
    
    fig, ax = plt.subplots(figsize=(12, 10))
    
    # Apply cuts and get filtered data
    cuts = plot_def.get('cuts')
    data_dict = apply_cuts(tree, cuts, plot_def['branches'])
    x_flat = data_dict[x_branch]
    y_flat = data_dict[y_branch]
    
    # Ensure same length (should be the case for supercluster data)
    min_len = min(len(x_flat), len(y_flat))
    x_flat = x_flat[:min_len]
    y_flat = y_flat[:min_len]
    
    # Create 2D histogram
    if 'binning' in plot_def and plot_def['binning'] is not None:
        x_bins, y_bins = plot_def['binning']
    else:
        x_bins = y_bins = 50  # Default 50x50

    counts, x_edges, y_edges = np.histogram2d(x_flat.to_numpy(), y_flat.to_numpy(), bins=[x_bins, y_bins])

    # Plot as image - transpose counts for correct orientation
    im = ax.imshow(counts.T, origin='lower', aspect='auto',
                   extent=[x_edges[0], x_edges[-1], y_edges[0], y_edges[-1]],
                   cmap='viridis', norm = "log" if logy_scale else None)
    
    # Add colorbar
    cbar = plt.colorbar(im, ax=ax)
    cbar.set_label('SuperClusters')
    
    # Set labels and title
    ax.set_xlabel(plot_def['xlabel'])
    ax.set_ylabel(plot_def['ylabel'])
    title = f"{plot_def['ylabel']} vs {plot_def['xlabel']}"
    # ax.set_title(title.replace('[', '(').replace(']', ')'), pad=40)
    
    # Add CMS label
    hep.cms.label("Simulation", data=False, ax=ax)
    
    # Add statistics
    entries = len(x_flat)
    stats_text = f'Entries: {entries}'
    if cuts is not None:
        stats_text += f'\n(with cuts applied)'
    
    ax.text(0.02, 0.98, stats_text, transform=ax.transAxes, fontsize=10,
            verticalalignment='top',
            bbox=dict(boxstyle="round,pad=0.3", facecolor="white", alpha=0.8))
    
    # Save the plot
    # cuts_suffix = "_cut" if cuts is not None else ""
    output_path = output_dir / f"{name}.png"
    plt.tight_layout()
    plt.savefig(output_path, dpi=300, bbox_inches='tight')
    plt.close()
    
    print(f"Saved 2D plot: {output_path}")

def main():
    parser = argparse.ArgumentParser(description='Plot MTD SuperCluster validation from TTree')
    parser.add_argument('input_file', nargs='?', default='mtd_supercluster_validation.root',
                       help='Input ROOT file (default: mtd_supercluster_validation.root)')
    parser.add_argument('--output-dir', '-o', default='/eos/home-n/npalmeri/www/MTD/SuperCluster/sim_tests',
                       help='Output directory for plots')
    parser.add_argument('--config', '-c', default='plot_config.py',
                       help='Plot configuration file (default: plot_config.py)')
    parser.add_argument('--tree-name', '-t', default='mtdSimSuperClusterValidation/MTDSuperClusters',
                       help='Name of TTree to read (default: mtdSimSuperClusterValidation/MTDSuperClusters)')
    
    args = parser.parse_args()
    
    # Check if input file exists
    if not os.path.exists(args.input_file):
        print(f"Error: Input file {args.input_file} not found!")
        return 1
    
    # Create output directory
    output_dir = Path(args.output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)
    
    print(f"Reading TTree from {args.input_file}...")
    print(f"Looking for tree: {args.tree_name}")

    # Load plot configuration (import PLOT_DEFINITIONS from there)
    if not os.path.exists(args.config):
        print(f"Error: Configuration file {args.config} not found!")
        return 1
    PLOT_DEFINITIONS = import_module(args.config.replace('.py', '').replace('/', '.')).PLOT_DEFINITIONS
    
    try:
        # Open ROOT file and get the tree
        with uproot.open(args.input_file) as file:
            # Check if tree exists
            if args.tree_name not in file:
                print(f"Error: Tree '{args.tree_name}' not found in ROOT file!")
                print("Available objects:")
                for key in file.keys():
                    print(f"  {key}")
                return 1
            
            tree = file[args.tree_name]
            
            # Print tree info
            print(f"Tree has {tree.num_entries} entries")
            print("Available branches:")
            for branch_name in tree.keys():
                print(f"  {branch_name}")
            
            # Check that required branches exist
            required_branches = set()
            for plot_def in PLOT_DEFINITIONS:
                required_branches.update(plot_def['branches'])
            
            missing_branches = required_branches - set(tree.keys())
            if missing_branches:
                print(f"\nError: Missing required branches: {missing_branches}")
                return 1
            
            print(f"\nCreating {len(PLOT_DEFINITIONS)} plots...")
            
            # Create plots for each definition
            for i, plot_def in enumerate(PLOT_DEFINITIONS):
                try:
                    print(f"\n[{i+1}/{len(PLOT_DEFINITIONS)}] Creating plot: {plot_def['name']}")
                    
                    if len(plot_def['branches']) == 1:
                        # 1D plot
                        plot_1d_from_tree(tree, plot_def, output_dir)
                    elif len(plot_def['branches']) == 2:
                        # 2D plot
                        plot_2d_from_tree(tree, plot_def, output_dir)
                    else:
                        print(f"Error: Plot {plot_def['name']} has {len(plot_def['branches'])} branches - only 1D and 2D plots supported")
                        continue
                        
                except Exception as e:
                    print(f"Error creating plot {plot_def['name']}: {e}")
                    continue
    
    except Exception as e:
        print(f"Error reading ROOT file: {e}")
        return 1
    
    print(f"\nAll plots saved to: {output_dir}")
    print("Done!")
    
    return 0

if __name__ == "__main__":
    exit(main())
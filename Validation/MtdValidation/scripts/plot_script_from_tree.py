#!/usr/bin/env python3
"""
Plot script for MTD MergedCluster validation from TTree
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
            if branch == 'simmc_n':
                result[branch] = tree[branch].array(library='np')
            else:
                data = tree[branch].array(library='ak')
                result[branch] = np.concatenate([event_data for event_data in data if len(event_data) > 0])
        return result
    
    # Read all needed branches including those used for cuts
    all_branches = set(branches_needed) | set(cuts.keys())
    branch_data = {}
    
    for branch in all_branches:
        if branch == 'simmc_n':
            branch_data[branch] = tree[branch].array(library='np')
        else:
            data = tree[branch].array(library='ak')
            branch_data[branch] = np.concatenate([event_data for event_data in data if len(event_data) > 0])
    
    # Apply cuts sequentially - ensure all arrays are numpy arrays
    first_mergedcluster_branch = None
    for branch in all_branches:
        if branch != 'simmc_n':
            first_mergedcluster_branch = branch
            break
    
    if first_mergedcluster_branch is None:
        # Only simmc_n branch, return as-is
        result = {}
        for branch in branches_needed:
            result[branch] = branch_data[branch]
        return result
    
    mask = np.ones(len(branch_data[first_mergedcluster_branch]), dtype=bool)
    
    for cut_branch, cut_func in cuts.items():
        if cut_branch in branch_data and cut_branch != 'simmc_n':
            # Convert to numpy array and apply cut
            cut_data = np.array(branch_data[cut_branch])
            cut_result = cut_func(cut_data)
            # Ensure cut_result is a numpy array
            if hasattr(cut_result, 'to_numpy'):
                cut_result = cut_result.to_numpy()
            else:
                cut_result = np.array(cut_result)
            mask = mask & cut_result
    
    # Apply mask to all branches except simmc_n (which is event-level)
    result = {}
    for branch in branches_needed:
        if branch == 'simmc_n':
            result[branch] = branch_data[branch]
        else:
            result[branch] = branch_data[branch][mask]
    
    return result

def plot_1d_from_tree(tree, plot_def, output_dir):
    """
    Create and plot a 1D histogram from one or more TTree plot definitions.

    plot_def may be a single plot definition (dict) or a list/tuple of plot definition dicts
    to be overlayed. Each definition should follow the existing format. When overlaying
    we compute a common binning (from the primary definition if provided, otherwise
    from the combined data) and draw each dataset with a different color and legend.
    """
    # Support passing either a single definition or a list of definitions to overlay
    if isinstance(plot_def, (list, tuple)):
        plot_defs = list(plot_def)
    else:
        plot_defs = [plot_def]

    # Use the first definition as the primary one for labels / default binning
    primary = plot_defs[0]
    name = primary['name']
    branch = primary['branches'][0]

    # Prepare figure
    fig, ax = plt.subplots(figsize=(10, 8))

    # Collect data for each definition separately (apply their own cuts)
    all_data = []
    all_type_info = []  # None or array of cluster types
    entries_list = []
    cuts_applied = False

    for pdef in plot_defs:
        # Only 1D plots supported for overlay
        if len(pdef['branches']) != 1:
            print(f"Error: Overlaid plot '{pdef['name']}' is not 1D. Skipping.")
            all_data.append(np.array([]))
            all_type_info.append(None)
            entries_list.append(0)
            continue

        pbranch = pdef['branches'][0]
        branches_needed = [pbranch]
        type_breakdown = pdef.get('type_breakdown', False)
        if type_breakdown and 'simmc_clusterType' not in branches_needed:
            branches_needed.append('simmc_clusterType')

        cuts = pdef.get('cuts')
        if cuts is not None:
            cuts_applied = True

        data_dict = apply_cuts(tree, cuts, branches_needed)

        # If the requested branch was not present in data_dict fall back to empty
        data = data_dict.get(pbranch, np.array([]))
        data = np.ravel(data)
        all_data.append(data)

        if type_breakdown:
            t = data_dict.get('simmc_clusterType')
            if t is not None:
                t = np.ravel(t)
                # Make sure it's aligned in length with data
                min_len = min(len(data), len(t))
                t = t[:min_len]
            all_type_info.append(t)
        else:
            all_type_info.append(None)

        entries_list.append(len(data))

    # Compute common binning: prefer explicit primary binning, otherwise derive from combined data
    if 'binning' in primary and primary['binning'] is not None:
        bins = primary['binning']
    else:
        # derive from combined non-empty data
        combined = np.concatenate([d for d in all_data if len(d) > 0]) if any(len(d) > 0 for d in all_data) else np.array([])
        if combined.size > 0:
            bins = np.histogram_bin_edges(combined, bins='auto')
        else:
            bins = 50

    # Prepare colors from matplotlib cycle
    color_cycle = plt.rcParams['axes.prop_cycle'].by_key().get('color', None)
    if color_cycle is None:
        # color_cycle = ['#5790fc', '#e42536', '#7a21dd', '#9c9ca1', '#2ca02c']
        color_cycle = ["#5790fc","#f89c20","#e42536","#964a8b","#9c9ca1","#7a21dd"]

    # Plot each dataset
    for idx, pdef in enumerate(plot_defs):
        data = all_data[idx]
        if data is None or len(data) == 0:
            # Still include legend entry with zero entries
            label = f"{pdef['name']} (0)"
            ax.plot([], [], label=label)
            continue

        # Ensure bins is a bin-edge array or integer accepted by numpy.histogram
        hist_bins = bins
        # If bins is an array of edges, use it directly; if it's an int, pass through
        counts, bin_edges = np.histogram(data, bins=hist_bins)
        bin_centers = (bin_edges[:-1] + bin_edges[1:]) / 2

        underflow = np.sum(data < bin_edges[0])
        overflow = np.sum(data > bin_edges[-1])
        uncertainties = np.sqrt(counts)

        base_color = color_cycle[idx % len(color_cycle)]
        # If this definition requests a type breakdown, draw the integrated (total)
        # histogram in black and rotate subtype colors from the color cycle.
        if all_type_info[idx] is not None:
            main_color = 'black'
        else:
            main_color = base_color

        label = f"{pdef['name']} ({len(data)})"

        # Plot error bars and step histogram (main / integrated)
        ax.errorbar(bin_centers, counts, yerr=uncertainties, fmt='o',
                    markersize=3, capsize=2, capthick=1,
                    linewidth=1.2, alpha=0.85, label=label, color=main_color)
        ax.hist(data, bins=hist_bins, histtype='step', linewidth=1.8,
                alpha=0.6, color=main_color, label='_nolegend_')

        # Optionally plot numbers above bins
        if pdef.get('plot_numbers', False):
            # small horizontal offset per overlay to avoid exact overlap
            n_overlays = len(plot_defs)
            offset_fraction = (idx - (n_overlays - 1) / 2) * 0.02  # adjust spacing
            # choose a base offset in data units (fraction of bin width)
            if np.isscalar(hist_bins):
                # integer bins: approximate width from bin_centers spacing
                bin_width = bin_centers[1] - bin_centers[0] if len(bin_centers) > 1 else 0.0
            else:
                # hist_bins is array of edges
                bin_width = bin_edges[1] - bin_edges[0] if len(bin_edges) > 1 else 0.0

            for xi, c in zip(bin_centers, counts):
                if c <= 0:
                    continue
                # compute horizontal position with a small fraction of bin width
                x_pos = xi + offset_fraction * (bin_width if bin_width != 0 else xi * 0.01)
                # compute vertical position slightly above the bar
                if primary.get('logy_scale', False):
                    # place labels at a small factor above the count on log scale
                    y_pos = max(c * 1.2, 0.1)
                else:
                    y_pos = c + max(0.02 * max(counts.max(), 1), 1)

                ax.text(x_pos, y_pos, f"{int(c)}",
                        ha='center', va='bottom', fontsize=15, rotation=45, color=main_color, alpha=0.9)

        # If this definition requests a type breakdown, draw them
        if all_type_info[idx] is not None:
            cluster_types = all_type_info[idx]
            # Align lengths
            min_len = min(len(data), len(cluster_types))
            data_for_types = data[:min_len]
            cluster_types = cluster_types[:min_len]

            cluster_type_names = {0: 'Primary', 1: 'Secondary', 2: 'Looper', 3: 'Calo backscatter'}

            unique_types = sorted(set(cluster_types))
            for sub_i, ctype in enumerate(unique_types):
                if ctype in cluster_type_names:
                    tmask = cluster_types == ctype
                    type_data = data_for_types[tmask]
                    if len(type_data) == 0:
                        continue
                    type_counts, _ = np.histogram(type_data, bins=hist_bins)
                    type_unc = np.sqrt(type_counts)
                    # rotate colors from the main color cycle so subtype colors vary
                    tcolor = color_cycle[(idx + sub_i) % len(color_cycle)]
                    tname = cluster_type_names[ctype]
                    ax.errorbar(bin_centers, type_counts, yerr=type_unc,
                                fmt='s', markersize=2, capsize=1, capthick=0.5,
                                linewidth=1, alpha=0.7, label=f"{pdef['name']}:{tname} ({len(type_data)})",
                                color=tcolor)
                    # print("PLOTTING DATA FOR ", pdef['name'], "type ", tname)
                    # print("MEAN AND STD OF DATA = ", np.mean(type_data), np.std(type_data))

                    ax.hist(type_data, bins=hist_bins, histtype='step', linewidth=1.2,
                            alpha=0.4, color=tcolor, label='_nolegend_')

    # Labels and scales from primary
    ax.set_xlabel(primary.get('xlabel', ''))
    ax.set_ylabel(primary.get('ylabel', ''))

    # Set log scale if requested on primary (keep simple)
    if primary.get('logy_scale', False):
        ax.set_yscale('log')
        # choose a safe bottom
        ax.set_ylim(bottom=0.1)
    if primary.get('logx_scale', False):
        ax.set_xscale('log')

    hep.cms.label(com=13.6, label="Preliminary", data=False, ax=ax)

    # Legend
    ax.legend(fontsize=9)

    # Statistics box: list entries for each overlayed dataset
    stats_lines = []
    for idx, pdef in enumerate(plot_defs):
        stats_lines.append(f"{pdef['name']}: {entries_list[idx]}")
    if cuts_applied:
        stats_lines.append('(with cuts applied)')
    stats_text = '\n'.join(stats_lines)
    ax.text(0.75, 0.85, stats_text, transform=ax.transAxes, fontsize=10,
            verticalalignment='top', bbox=dict(boxstyle="round,pad=0.3", facecolor="white", alpha=0.8))

    # Save the plot
    if len(plot_defs) > 1:
        overlap_suffix = '_overlap_' + '_'.join([p['name'] for p in plot_defs[1:]])
        output_name = f"{name}{overlap_suffix}.png"
    else:
        output_name = f"{name}.png"

    output_path = output_dir / output_name
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
    
    # Ensure same length (should be the case for mergedcluster data)
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
    cbar.set_label('MergedClusters')
    
    # Set labels and title
    ax.set_xlabel(plot_def['xlabel'])
    ax.set_ylabel(plot_def['ylabel'])
    title = f"{plot_def['ylabel']} vs {plot_def['xlabel']}"
    # ax.set_title(title.replace('[', '(').replace(']', ')'), pad=40)
    
    # Add CMS label
    hep.cms.label(com=13.6, label="Preliminary", data=False, ax=ax)
    
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
    parser = argparse.ArgumentParser(description='Plot MTD MergedCluster validation from TTree')
    parser.add_argument('input_file', nargs='?', default='mtd_mergedcluster_validation.root',
                       help='Input ROOT file (default: mtd_mergedcluster_validation.root)')
    parser.add_argument('--output-dir', '-o', default='/eos/home-n/npalmeri/www/MTD/MergedCluster/sim_tests',
                       help='Output directory for plots')
    parser.add_argument('--config', '-c', default='plot_config.py',
                       help='Plot configuration file (default: plot_config.py)')
    parser.add_argument('--tree-name', '-t', default='mtdMergedClusterValidation/MTDMergedClusters',
                       help='Name of TTree to read (default: mtdMergedClusterValidation/MTDMergedClusters)')
    
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
            plotted = set()
            name_to_def = {pd['name']: pd for pd in PLOT_DEFINITIONS}

            inline_counter = 0
            for i, plot_def in enumerate(PLOT_DEFINITIONS):
                pname = plot_def['name']
                if pname in plotted:
                    continue

                try:
                    print(f"\n[{i+1}/{len(PLOT_DEFINITIONS)}] Creating plot: {pname}")

                    # Check for overlap specification: can be a list of names or inline defs (dicts)
                    overlap_list = plot_def.get('overlap')
                    if overlap_list:
                        # Gather plot defs to overlay (primary + overlaps)
                        overlay_defs = [plot_def]
                        for oname in overlap_list:
                            # support strings referencing named definitions
                            if isinstance(oname, str):
                                if oname in name_to_def:
                                    overlay_defs.append(name_to_def[oname])
                                    # mark the named overlay as plotted to avoid duplicate standalone plotting
                                    plotted.add(oname)
                                else:
                                    print(f"Warning: Overlap target '{oname}' not found in PLOT_DEFINITIONS")
                            # support inline dict definitions so users can provide different cuts without creating a named def
                            elif isinstance(oname, dict):
                                # ensure it looks like a plot definition
                                if 'branches' not in oname:
                                    print(f"Warning: Inline overlap definition is missing 'branches' and will be skipped: {oname}")
                                    continue
                                # if no name provided, generate a unique inline name
                                if 'name' not in oname:
                                    oname = dict(oname)  # copy to avoid mutating user's object
                                    oname['name'] = f"inline_overlay_{inline_counter}_{oname['branches'][0]}"
                                    inline_counter += 1
                                overlay_defs.append(oname)
                            else:
                                print(f"Warning: Overlap entry must be a string name or a dict definition. Skipping: {oname}")

                        # Ensure all overlay defs are 1D
                        if any(len(pd['branches']) != 1 for pd in overlay_defs):
                            print(f"Error: One of the overlaid plots for '{pname}' is not 1D. Skipping overlay.")
                        else:
                            plot_1d_from_tree(tree, overlay_defs, output_dir)
                            plotted.add(pname)

                    else:
                        # No overlap; plot normally
                        if len(plot_def['branches']) == 1:
                            plot_1d_from_tree(tree, plot_def, output_dir)
                            plotted.add(pname)
                        elif len(plot_def['branches']) == 2:
                            plot_2d_from_tree(tree, plot_def, output_dir)
                            plotted.add(pname)
                        else:
                            print(f"Error: Plot {pname} has {len(plot_def['branches'])} branches - only 1D and 2D plots supported")
                            continue

                except Exception as e:
                    print(f"Error creating plot {pname}: {e}")
                    continue
    
    except Exception as e:
        print(f"Error reading ROOT file: {e}")
        return 1
    
    print(f"\nAll plots saved to: {output_dir}")
    print("Done!")
    
    return 0

if __name__ == "__main__":
    exit(main())
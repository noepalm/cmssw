#!/usr/bin/env python3
"""
Plot script for MTD MergedCluster validation from TTree
Creates maps of cluster energy vs iphi/ieta for events where topological clustering merges clusters
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

def get_mtd_ranges():
    """
    Return fixed MTD iPhi/iEta ranges
    """
    # # Fixed MTD geometry: iPhi 0-107 (108 total), iEta 0-95 (96 total)
    # return 0, 107, 0, 95
    return 1, 108, 1, 96

def create_cluster_map_combined(iphi_values, ieta_values, energy_values, time_values, clustertype_values, title_base, output_path, mtd_ranges, primary_energy=None, primary_eta=None, primary_phi=None, primary_pt=None):
    """
    Create a combined 2D histogram map showing energy, time, and cluster type weighted cluster positions side-by-side
    """
    if len(iphi_values) == 0:
        print(f"Warning: No data for {title_base}")
        return
    
    # Create larger figure with three subplots side by side
    fig, (ax1, ax2, ax3) = plt.subplots(1, 3, figsize=(40, 12))
    
    # Use fixed MTD ranges
    iphi_min, iphi_max, ieta_min, ieta_max = mtd_ranges
    
    # Create histogram bins
    iphi_bins = np.arange(iphi_min, iphi_max + 2)
    ieta_bins = np.arange(ieta_min, ieta_max + 2)
    
    # Energy-weighted histogram (left subplot)
    hist_energy, iphi_edges, ieta_edges = np.histogram2d(
        iphi_values, ieta_values, bins=[iphi_bins, ieta_bins], weights=energy_values
    )
    energy_min = np.min(energy_values) * 0.1 if len(energy_values) > 0 else 0
    hist_energy = np.where(hist_energy == 0, np.nan, hist_energy)
    
    im1 = ax1.imshow(hist_energy.T, origin='lower', aspect='auto',
                     extent=[iphi_edges[0], iphi_edges[-1], ieta_edges[0], ieta_edges[-1]],
                     cmap='viridis', interpolation='nearest', vmin=energy_min)
    
    # Time-weighted histogram (middle subplot)
    hist_time, _, _ = np.histogram2d(
        iphi_values, ieta_values, bins=[iphi_bins, ieta_bins], weights=time_values
    )
    time_min = np.min(time_values) * 0.1 if len(time_values) > 0 else 0
    hist_time = np.where(hist_time == 0, np.nan, hist_time)
    
    im2 = ax2.imshow(hist_time.T, origin='lower', aspect='auto',
                     extent=[iphi_edges[0], iphi_edges[-1], ieta_edges[0], ieta_edges[-1]],
                     cmap='plasma', interpolation='nearest', vmin=time_min)
    
    # print("\tDEBUG: Cluster types: ")
    # for idx in range(len(clustertype_values)):
    #     print(f"\t  Cluster {idx}: iphi = {iphi_values[idx]}, ieta = {ieta_values[idx]}, type={clustertype_values[idx]}")

    # print("Is there at least one primary cluster? (Cluster type 0)", any(ct == 0 for ct in clustertype_values))

    # Cluster type histogram (right subplot)
    hist_type, _, _ = np.histogram2d(
        iphi_values, ieta_values, bins=[iphi_bins, ieta_bins], weights=np.array(clustertype_values) + 1
    )

    # For cluster type, we want to show the actual values, so use discrete colormap
    hist_type = np.where(hist_type == 0, np.nan, hist_type)
    
    # Create a discrete colormap for cluster types
    from matplotlib.colors import ListedColormap, BoundaryNorm
    cluster_colors = ['#e42536', '#5790fc', '#7a21dd', '#9c9ca1']  # Red, Blue, Purple, Gray
    cluster_cmap = ListedColormap(cluster_colors)
    cluster_norm = BoundaryNorm(boundaries=[0.5, 1.5, 2.5, 3.5, 4.5], ncolors=4)
    
    im3 = ax3.imshow(hist_type.T, origin='lower', aspect='auto',
                     extent=[iphi_edges[0], iphi_edges[-1], ieta_edges[0], ieta_edges[-1]],
                     cmap=cluster_cmap, norm=cluster_norm, interpolation='nearest',
                     alpha=1)
    
    # Plot primary clusters (type 0) on top with full opacity for visibility
    primary_indices = [i for i, ct in enumerate(clustertype_values) if ct == 0]
    if primary_indices:
        primary_iphi = [iphi_values[i] for i in primary_indices]
        primary_ieta = [ieta_values[i] for i in primary_indices]
        primary_weights = [1 for _ in primary_indices]  # Weight of 1 for primary clusters (shifted becomes 1)
        
        hist_primary, _, _ = np.histogram2d(
            primary_iphi, primary_ieta, bins=[iphi_bins, ieta_bins], weights=primary_weights
        )
        hist_primary = np.where(hist_primary == 0, np.nan, hist_primary)
        
        # Plot primary clusters on top with full opacity
        ax3.imshow(hist_primary.T, origin='lower', aspect='auto',
                   extent=[iphi_edges[0], iphi_edges[-1], ieta_edges[0], ieta_edges[-1]],
                   cmap=ListedColormap(['#e42536']), interpolation='nearest',
                   alpha=1.0, vmin=0.5, vmax=1.5)
    
    # Configure all subplots
    subplot_configs = [
        (ax1, im1, "Energy", "Energy [MeV]"),
        (ax2, im2, "Time", "Time [ns]"),
        (ax3, im3, "Cluster Type", "Cluster Type")
    ]
    
    for ax, im, title_suffix, colorbar_label in subplot_configs:
        # Add colorbar
        cbar = plt.colorbar(im, ax=ax)
        cbar.set_label(colorbar_label, fontsize=12)
        
        # Special handling for cluster type colorbar
        if title_suffix == "Cluster Type":
            cbar.set_ticks([1, 2, 3, 4])
            cbar.set_ticklabels(['Primary', 'Secondary', 'Looper', 'Backscatter'])
        
        # Set labels and title
        ax.set_xlabel('iPhi', fontsize=12)
        ax.set_ylabel('iEta', fontsize=12)
        
        # Add primary photon information to title if available
        enhanced_title = f"{title_base} ({title_suffix})"
        if primary_energy is not None and primary_energy > 0:
            enhanced_title += f"\nPhoton: E={primary_energy:.1f} GeV, pT={primary_pt:.1f} GeV, eta={primary_eta:.2f}, phi={primary_phi:.2f}"
        
        ax.set_title(enhanced_title, fontsize=14, pad=20)
        
        # Add grid with better spacing for readability
        ax.set_xticks(np.arange(iphi_min, iphi_max + 1))
        ax.set_yticks(np.arange(ieta_min, ieta_max + 1))
        
        # Make tick labels smaller and give them more space
        ax.tick_params(axis='x', labelsize=9, rotation=60, pad=8)
        ax.tick_params(axis='y', labelsize=10, rotation=0, pad=5)
        
        ax.grid(True, alpha=0.3)
        
        # Add statistics
        n_clusters = len(iphi_values)
        total_energy = sum(energy_values)
        stats_text = f'Clusters: {n_clusters}\nTotal Energy: {total_energy:.2f} MeV'
        ax.text(0.02, 0.98, stats_text, transform=ax.transAxes, fontsize=10,
                verticalalignment='top',
                bbox=dict(boxstyle="round,pad=0.3", facecolor="white", alpha=0.8))
    
    # Save the plot with more spacing
    plt.tight_layout()
    plt.savefig(output_path, dpi=300, bbox_inches='tight')
    plt.close()
    
    print(f"Saved combined map: {output_path}")

def process_event_maps(history_tree, topo_tree, event_idx, output_dir, mtd_ranges):
    """
    Process a single event and create cluster maps
    """
    # Get data for this specific event - using new per-cluster branch names
    evt_history = {
        'evt_number': history_tree['evt_number'].array(library='np')[event_idx],
        'sc_n': history_tree['sc_n'].array(library='np')[event_idx],
        'sc_energy': history_tree['sc_energy'].array(library='ak')[event_idx],
        'sc_nClusters': history_tree['sc_nClusters'].array(library='ak')[event_idx],
        'sc_iphi_perCluster': history_tree['sc_iphi_perCluster'].array(library='ak')[event_idx],
        'sc_ieta_perCluster': history_tree['sc_ieta_perCluster'].array(library='ak')[event_idx],
        'sc_energy_perCluster': history_tree['sc_energy_perCluster'].array(library='ak')[event_idx],
        'sc_time_perCluster': history_tree['sc_time_perCluster'].array(library='ak')[event_idx],
        'sc_clusterType': history_tree['sc_clusterType'].array(library='ak')[event_idx],
        'sc_primary_energy': history_tree['sc_primary_energy'].array(library='ak')[event_idx],
        'sc_primary_et': history_tree['sc_primary_et'].array(library='ak')[event_idx],
        'sc_primary_phi': history_tree['sc_primary_phi'].array(library='ak')[event_idx],
        'sc_primary_eta': history_tree['sc_primary_eta'].array(library='ak')[event_idx]
    }
    
    evt_topo = {
        'evt_number': topo_tree['evt_number'].array(library='np')[event_idx],
        'sc_n': topo_tree['sc_n'].array(library='np')[event_idx],
        'sc_energy': topo_tree['sc_energy'].array(library='ak')[event_idx],
        'sc_nClusters': topo_tree['sc_nClusters'].array(library='ak')[event_idx],
        'sc_iphi_perCluster': topo_tree['sc_iphi_perCluster'].array(library='ak')[event_idx],
        'sc_ieta_perCluster': topo_tree['sc_ieta_perCluster'].array(library='ak')[event_idx],
        'sc_energy_perCluster': topo_tree['sc_energy_perCluster'].array(library='ak')[event_idx],
        'sc_time_perCluster': topo_tree['sc_time_perCluster'].array(library='ak')[event_idx],
        'sc_clusterType': topo_tree['sc_clusterType'].array(library='ak')[event_idx],
        'sc_primary_energy': topo_tree['sc_primary_energy'].array(library='ak')[event_idx],
        'sc_primary_et': topo_tree['sc_primary_et'].array(library='ak')[event_idx],
        'sc_primary_phi': topo_tree['sc_primary_phi'].array(library='ak')[event_idx],
        'sc_primary_eta': topo_tree['sc_primary_eta'].array(library='ak')[event_idx]
    }
    
    # Create event directory
    event_dir = output_dir / str(event_idx)
    event_dir.mkdir(exist_ok=True)
    
    # Get the actual event number for display
    event_number = evt_history['evt_number']
    
    print(f"Processing event {event_idx} (event number {event_number}):")
    print(f"  History-only: {evt_history['sc_n']} mergedclusters")
    print(f"  Topo+History: {evt_topo['sc_n']} mergedclusters")
    
    # Debug prints for topological+historical clustering mergedclusters with >=2 clusters
    topo_merged_info = []
    for sc_idx in range(len(evt_topo['sc_energy'])):
        n_clusters = evt_topo['sc_nClusters'][sc_idx]
        if n_clusters >= 2:
            # Now we have per-cluster data - show all cluster coordinates
            iphi_list = list(evt_topo['sc_iphi_perCluster'][sc_idx])
            ieta_list = list(evt_topo['sc_ieta_perCluster'][sc_idx])
            energy_list = list(evt_topo['sc_energy_perCluster'][sc_idx])
            time_list = list(evt_topo['sc_time_perCluster'][sc_idx])
            clustertype_list = list(evt_topo['sc_clusterType'][sc_idx])
            total_energy = evt_topo['sc_energy'][sc_idx]
            topo_merged_info.append({
                'sc_idx': sc_idx,
                'n_clusters': n_clusters,
                'iphi_list': iphi_list,
                'ieta_list': ieta_list,
                'energy_list': energy_list,
                'time_list': time_list,
                'clustertype_list': clustertype_list,
                'total_energy': total_energy
            })
    
    print(f"  Topo+History merged mergedclusters ({len(topo_merged_info)}):")
    for info in topo_merged_info:
        print(f"    SC {info['sc_idx']}: {info['n_clusters']} clusters, Total E={info['total_energy']:.2f} MeV")
        cluster_type_names = {0: 'Primary', 1: 'Secondary', 2: 'Looper', 3: 'Calo backscatter'}
        for i in range(len(info['iphi_list'])):
            cluster_type_name = cluster_type_names.get(info['clustertype_list'][i], f"Unknown({info['clustertype_list'][i]})")
            print(f"      Cluster {i}: iPhi={info['iphi_list'][i]}, iEta={info['ieta_list'][i]}, E={info['energy_list'][i]:.2f} MeV, t={info['time_list'][i]:.2f} ns, type={cluster_type_name}")
    
    if len(topo_merged_info) == 0:
        print("    No merged mergedclusters found (this shouldn't happen)")
        return False
    
    # Process history-only clustering - plot each mergedcluster with individual cluster positions
    for sc_idx in range(len(evt_history['sc_energy'])):
        n_clusters = evt_history['sc_nClusters'][sc_idx]
        if n_clusters <= 0:
            continue
            
        # Now we can plot individual clusters within the mergedcluster
        iphi_vals = list(evt_history['sc_iphi_perCluster'][sc_idx])
        ieta_vals = list(evt_history['sc_ieta_perCluster'][sc_idx])
        energy_vals = list(evt_history['sc_energy_perCluster'][sc_idx])
        time_vals = list(evt_history['sc_time_perCluster'][sc_idx])
        clustertype_vals = list(evt_history['sc_clusterType'][sc_idx])
        
        # Create combined energy, time, and cluster type map
        title_base = f"Event {event_idx} (#{event_number}) - History Only - SC {sc_idx}\n({n_clusters} clusters, {evt_history['sc_energy'][sc_idx]:.2f} MeV total)"
        output_path = event_dir / f"history_sc_{sc_idx}.png"
        
        # Get primary particle info for this mergedcluster
        primary_energy = evt_history['sc_primary_energy'][sc_idx] if evt_history['sc_primary_energy'][sc_idx] > 0 else None
        primary_eta = evt_history['sc_primary_eta'][sc_idx] if primary_energy is not None else None
        primary_phi = evt_history['sc_primary_phi'][sc_idx] if primary_energy is not None else None
        primary_pt = evt_history['sc_primary_et'][sc_idx] if primary_energy is not None else None
        
        create_cluster_map_combined(iphi_vals, ieta_vals, energy_vals, time_vals, clustertype_vals, title_base, output_path, mtd_ranges, primary_energy, primary_eta, primary_phi, primary_pt)
    
    # Process topo+history clustering - only mergedclusters with >=2 clusters, showing individual cluster positions
    merged_sc_count = 0
    for sc_idx in range(len(evt_topo['sc_energy'])):
        n_clusters = evt_topo['sc_nClusters'][sc_idx]
        if n_clusters < 2:
            continue
            
        merged_sc_count += 1
        
        # Plot individual clusters within the merged mergedcluster
        iphi_vals = list(evt_topo['sc_iphi_perCluster'][sc_idx])
        ieta_vals = list(evt_topo['sc_ieta_perCluster'][sc_idx])
        energy_vals = list(evt_topo['sc_energy_perCluster'][sc_idx])
        time_vals = list(evt_topo['sc_time_perCluster'][sc_idx])
        clustertype_vals = list(evt_topo['sc_clusterType'][sc_idx])
        
        # Create combined energy, time, and cluster type map
        title_base = f"Event {event_idx} (#{event_number}) - Topo+History - SC {sc_idx}\n({n_clusters} clusters, {evt_topo['sc_energy'][sc_idx]:.2f} MeV total)"
        output_path = event_dir / f"topo_sc_{sc_idx}_merged.png"
        
        # Get primary particle info for this mergedcluster
        primary_energy = evt_topo['sc_primary_energy'][sc_idx] if evt_topo['sc_primary_energy'][sc_idx] > 0 else None
        primary_eta = evt_topo['sc_primary_eta'][sc_idx] if primary_energy is not None else None
        primary_phi = evt_topo['sc_primary_phi'][sc_idx] if primary_energy is not None else None
        primary_pt = evt_topo['sc_primary_et'][sc_idx] if primary_energy is not None else None
        
        create_cluster_map_combined(iphi_vals, ieta_vals, energy_vals, time_vals, clustertype_vals, title_base, output_path, mtd_ranges, primary_energy, primary_eta, primary_phi, primary_pt)
    
    print(f"  Created maps for {merged_sc_count} merged mergedclusters")
    return merged_sc_count > 0

def main():
    parser = argparse.ArgumentParser(description='Create cluster maps for events with topological merging')
    parser.add_argument('--history-file', default='mtd_mergedcluster_validation_history.root',
                       help='History-only clustering ROOT file')
    parser.add_argument('--topo-file', default='mtd_mergedcluster_validation_topo+history.root',
                       help='Topological+history clustering ROOT file')
    parser.add_argument('--output-dir', '-o', default='/eos/home-n/npalmeri/www/MTD/MergedCluster/sim_tests',
                       help='Output directory for plots')
    parser.add_argument('--tree-name', '-t', default='mtdSimMergedClusterValidation/MTDMergedClusters',
                       help='Name of TTree to read')
    parser.add_argument('--max-events', default=10, type=int,
                       help='Maximum event number to process (default: 10)')
    parser.add_argument('--single-event', type=int, default=None,
                       help='Process only a single event by index (overrides max-events)')
    
    args = parser.parse_args()
    
    # Check if input files exist
    if not os.path.exists(args.history_file):
        print(f"Error: History file {args.history_file} not found!")
        return 1
    
    if not os.path.exists(args.topo_file):
        print(f"Error: Topology file {args.topo_file} not found!")
        return 1
    
    # Create output directory
    output_dir = Path(args.output_dir) / "event"
    output_dir.mkdir(parents=True, exist_ok=True)
    
    print(f"Reading trees...")
    print(f"History file: {args.history_file}")
    print(f"Topology file: {args.topo_file}")
    print(f"Output directory: {output_dir}")
    
    try:
        # Open both ROOT files
        with uproot.open(args.history_file) as history_file, \
             uproot.open(args.topo_file) as topo_file:
            
            # Get trees
            if args.tree_name not in history_file:
                print(f"Error: Tree '{args.tree_name}' not found in history file!")
                return 1
            
            if args.tree_name not in topo_file:
                print(f"Error: Tree '{args.tree_name}' not found in topology file!")
                return 1
            
            history_tree = history_file[args.tree_name]
            topo_tree = topo_file[args.tree_name]
            
            # Check that trees have same number of entries
            n_entries_history = history_tree.num_entries
            n_entries_topo = topo_tree.num_entries
            
            if n_entries_history != n_entries_topo:
                print(f"Warning: Trees have different number of entries!")
                print(f"  History: {n_entries_history}")
                print(f"  Topology: {n_entries_topo}")
            
            n_entries = min(n_entries_history, n_entries_topo)
            
            # Handle single event processing
            if args.single_event is not None:
                if args.single_event >= n_entries:
                    print(f"Error: Single event index {args.single_event} is out of range (0-{n_entries-1})")
                    return 1
                print(f"Processing single event at index {args.single_event}")
                events_to_check = [args.single_event]
            else:
                n_entries = min(n_entries, args.max_events)
                print(f"Processing {n_entries} events...")
                events_to_check = list(range(n_entries))
            
            # Check required branches exist - updated for new per-cluster branch names
            required_branches = ['evt_number', 'sc_n', 'sc_energy', 'sc_nClusters', 'sc_iphi_perCluster', 'sc_ieta_perCluster', 'sc_energy_perCluster', 'sc_time_perCluster', 'sc_clusterType', 'sc_primary_energy', 'sc_primary_et', 'sc_primary_phi', 'sc_primary_eta']
            for branch in required_branches:
                if branch not in history_tree.keys():
                    print(f"Error: Branch '{branch}' not found in history tree!")
                    return 1
                if branch not in topo_tree.keys():
                    print(f"Error: Branch '{branch}' not found in topology tree!")
                    return 1
            
            # Get fixed MTD ranges
            mtd_ranges = get_mtd_ranges()
            print(f"Using fixed MTD ranges: iPhi [1, 108], iEta [1, 96]")
            
            # Read topology clustering data to find interesting events
            if args.single_event is not None:
                # For single event, read only that event's data
                topo_nClusters = topo_tree['sc_nClusters'].array(library='ak', entry_start=args.single_event, entry_stop=args.single_event+1)
            else:
                # For multiple events, read up to max_events
                max_entries = min(n_entries_history, n_entries_topo, args.max_events)
                topo_nClusters = topo_tree['sc_nClusters'].array(library='ak')[:max_entries]
            
            events_with_merging = []
            for i, event_idx in enumerate(events_to_check):
                # Check if this event has any mergedclusters with >=2 clusters in topology
                event_clusters = topo_nClusters[i] if args.single_event is not None else topo_nClusters[event_idx]
                has_merged = any(n_clus >= 2 for n_clus in event_clusters)
                
                if has_merged:
                    events_with_merging.append(event_idx)
            
            if args.single_event is not None:
                if len(events_with_merging) == 0:
                    print(f"Event {args.single_event} does not have topological merging. Processing anyway.")
                    events_to_process = [args.single_event]
                else:
                    print(f"Event {args.single_event} has topological merging. Processing.")
                    events_to_process = events_with_merging
            else:
                print(f"Found {len(events_with_merging)} total events with topological merging")
                
                if len(events_with_merging) == 0:
                    print("No events with merging found. Exiting.")
                    return 0
                
                # Limit events to max_events (by event index, not count of merged events)
                events_to_process = [idx for idx in events_with_merging if idx < args.max_events]
                print(f"Will process {len(events_to_process)} events with merging in first {args.max_events} events")
            
            # Process each interesting event
            processed_count = 0
            for event_idx in events_to_process:
                try:
                    success = process_event_maps(history_tree, topo_tree, event_idx, output_dir, mtd_ranges)
                    if success:
                        processed_count += 1
                except Exception as e:
                    print(f"Error processing event {event_idx}: {e}")
                    continue
            
            print(f"\nSuccessfully processed {processed_count} events with merging")
            print(f"Maps saved to: {output_dir}")
    
    except Exception as e:
        print(f"Error: {e}")
        return 1
    
    return 0

if __name__ == "__main__":
    exit(main())


import uproot
import matplotlib.pyplot as plt
import mplhep as hep
import numpy as np
import awkward as ak
import os
import argparse
from pathlib import Path
from importlib import import_module

def main():
    # Set matplotlib style
    plt.style.use(hep.style.CMS)

    parser = argparse.ArgumentParser(description='Plot complicated MTD MergedCluster plots')
    parser.add_argument('input_file', nargs='?', default='../test/tree_vali_dev.root',
                        help='Input ROOT file (default: ../test/tree_vali_dev.root)')
    parser.add_argument('--output-dir', '-o', default='/eos/home-n/npalmeri/www/MTD/MergedCluster/forPresentation_251024/photonGun',
                        help='Output directory for plots')

    args = parser.parse_args()

    # Check if input file exists
    if not os.path.exists(args.input_file):
        print(f"Error: Input file {args.input_file} not found!")
        return

    # Create output directory
    output_dir = Path(args.output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)

    print(f"Reading data from {args.input_file}...")

    f = uproot.open(args.input_file)
    tree = f["mtdMergedClusterValidation/MTDMergedClusters"]
    arrays = tree.arrays()

    # only select merged clusters with 2 clusters inside
    nclusters_mask = arrays["simmc_nClusters"] == 2

    # take ieta, iphi
    ieta = arrays["simmc_ieta_perCluster"][nclusters_mask]
    iphi = arrays["simmc_iphi_perCluster"][nclusters_mask]

    # only keep entries with exactly 2 clusters
    # (NOTE: also selecting only first cluster per event for simplicity)
    ieta = ieta[ak.sum(nclusters_mask, axis = 1) > 0]
    iphi = iphi[ak.sum(nclusters_mask, axis = 1) > 0]

    # convert arrays to SIGNED ints
    ieta = ak.values_astype(ieta, np.int16)
    iphi = ak.values_astype(iphi, np.int16)

    # MANUALLY COMPUTE delta phi, delta eta
    # iterate over all events with 2 clusters
    dphi = []
    deta = []

    for e_eta, e_phi in zip(ieta, iphi):
        for cluster_eta, cluster_phi in zip(e_eta, e_phi):
            dphi.append(cluster_phi[0] - cluster_phi[1])
            deta.append(cluster_eta[0] - cluster_eta[1])

    # plot 2D histogram of dphi vs deta using integer-aligned, larger bins
    # convert to numpy integer arrays
    deta = np.array(deta, dtype=int)
    dphi = np.array(dphi, dtype=int)

    # determine bin edges so each integer value gets its own bin (width=1)
    min_deta, max_deta = -10, 10
    min_dphi, max_dphi = -10, 10
    # create bin edges that center on integers and span the full observed range
    x_bins = np.arange(min_deta - 0.5, max_deta + 1.5, 1.0)
    y_bins = np.arange(min_dphi - 0.5, max_dphi + 1.5, 1.0)

    plt.figure(figsize=(10, 8.5))
    plt.hist2d(dphi, deta, bins=[x_bins, y_bins])
    plt.colorbar(label='Counts')
    # enable latex labels
    # plt.rcParams["text.usetex"] = True
    plt.xlabel(r'$\Delta i \phi$')
    plt.ylabel(r'$\Delta i \eta$')

    plt.grid(True, linestyle='--', alpha=0.3, color="white")

    # write number of entries on the plot
    n_entries = len(deta)
    plt.text(0.95, 0.95, f'Entries: {n_entries}', horizontalalignment='right', verticalalignment='top', transform=plt.gca().transAxes, color='white', fontsize=12, bbox=dict(facecolor='black', alpha=0.5, pad=5))

    # change ticks to be at integers only
    plt.xticks(np.arange(min_dphi, max_dphi + 1, 1))
    plt.yticks(np.arange(min_deta, max_deta + 1, 1))

    # add cms label
    hep.cms.label("Preliminary", data=False, com = 13.6)

    # save figure
    output_path = output_dir / "deta_dphi_2clusters.png"
    plt.tight_layout()
    plt.savefig(output_path)
    print(f"Saved plot: {output_path}")    # 2. PLOT CLUSTER TYPES for those

    types = arrays["simmc_clusterType"][nclusters_mask]    # only select merged clusters with 2 clusters inside
    types = types[ak.sum(nclusters_mask, axis = 1) > 0]

    # plot cluster type 1 vs cluster type 2
    type1 = []
    type2 = []

    for e_types in types:
        for cluster_type in e_types:
            type1.append(cluster_type[0])
            type2.append(cluster_type[1])

    plt.figure(figsize=(10, 8.5))
    plt.hist2d(type1, type2, bins=[np.arange(-0.5, 4.5, 1), np.arange(-0.5, 4.5, 1)])
    plt.colorbar(label='Counts')
    plt.xlabel('Cluster Type 1')
    plt.ylabel('Cluster Type 2')

    plt.grid(True, linestyle='--', alpha=0.3, color="white")
    # write number of entries on the plot
    n_entries = len(type1)
    plt.text(0.95, 0.95, f'Entries: {n_entries}', horizontalalignment='right', verticalalignment='top', transform=plt.gca().transAxes,
                color='white', fontsize=12, bbox=dict(facecolor='black', alpha=0.5, pad=5))
    # set ticks to be at integers only
    plt.xticks(np.arange(0, 4, 1))
    plt.yticks(np.arange(0, 4, 1))

    # write number of entries in each BIN
    for i in range(4):
        for j in range(4):
            # count number of entries in each bin
            count = np.sum((np.array(type1) == i) & (np.array(type2) == j))
            if count > 0:
                plt.text(i, j, str(count), color='white', fontsize=15, ha='center', va='center')


    # add cms label
    hep.cms.label("Preliminary", data=False, com = 13.6)
    # save figure
    output_path = output_dir / "clustertype_2clusters.png"
    plt.tight_layout()
    plt.savefig(output_path)
    print(f"Saved plot: {output_path}")

    ### plot the same as above, but with fractions instead of absolute numbers

    plt.figure(figsize=(10, 8.5))
    plt.hist2d(type1, type2, bins=[np.arange(-0.5, 4.5, 1), np.arange(-0.5, 4.5, 1)])
    plt.colorbar(label='Counts')
    plt.xlabel('Cluster Type 1')
    plt.ylabel('Cluster Type 2')

    plt.grid(True, linestyle='--', alpha=0.3, color="white")
    # write number of entries on the plot
    n_entries = len(type1)
    plt.text(0.95, 0.95, f'Entries: {n_entries}', horizontalalignment='right', verticalalignment='top', transform=plt.gca().transAxes,
                color='white', fontsize=12, bbox=dict(facecolor='black', alpha=0.5, pad=5))
    # set ticks to be at integers only
    plt.xticks(np.arange(0, 4, 1))
    plt.yticks(np.arange(0, 4, 1))

    # write number of entries in each BIN
    for i in range(4):
        for j in range(4):
            # count number of entries in each bin
            count = np.sum((np.array(type1) == i) & (np.array(type2) == j)) / n_entries * 100
            if count > 0:
                plt.text(i, j, f"{count:.1f}%", color='white', fontsize=20, ha='center', va='center')

    # add cms label
    hep.cms.label("Preliminary", data=False, com = 13.6)
    # save figure
    output_path = output_dir / "clustertype_fractions_2clusters.png"
    plt.tight_layout()
    plt.savefig(output_path)
    print(f"Saved plot: {output_path}")

    ### 3. plot time, energy by cluster type for MERGED CLUSTERS ONLY
    clusterTypes = arrays["simmc_clusterType"]
    nClusters = arrays["simmc_nClusters"]
    energy = arrays["simmc_energy"]
    time = arrays["simmc_time"]

    mask = nClusters > 1
    good_evts = ak.num(energy[mask]) > 0

    energy_to_plot = energy[mask][good_evts]
    time_to_plot = time[mask][good_evts]
    clusterType_to_plot = ak.firsts(clusterTypes[mask][good_evts], axis = 2)

    # Define cluster type labels
    cluster_labels = {0: 'primary', 1: 'secondary', 2: 'loopers', 3: 'backscatter'}

    # plot energy by cluster type
    plt.figure(figsize=(10, 8.5))
    energy_bins = np.logspace(-3, 3, 15)
    for ctype in range(4):
        sel = clusterType_to_plot == ctype
        data = ak.flatten(energy_to_plot[sel])
        counts, bin_edges = np.histogram(data, bins=energy_bins)
        # Use geometric mean for log-scale bin centers
        bin_centers = np.sqrt(bin_edges[:-1] * bin_edges[1:])
        errors = np.sqrt(counts)
        # Plot histogram and get the color
        hist_plot = plt.hist(data, bins=energy_bins, histtype='step', label=cluster_labels[ctype], linewidth=1.5, alpha=0.8)
        hist_color = hist_plot[2][0].get_edgecolor()
        # Plot error bars with matching color and central marker
        plt.errorbar(bin_centers, counts, yerr=errors, fmt='o', markersize=3, elinewidth=1.2, capsize=2.5, capthick=1.2, color=hist_color, alpha=0.9)
    plt.xscale('log')
    plt.xlabel('MergedCluster Energy [MeV]')
    plt.ylabel('Entries')
    plt.legend()

    # add cms label
    hep.cms.label("Preliminary", data=False, com = 13.6)
    # save figure
    output_path = output_dir / "energy_byClusterType_reallyMergedOnly.png"
    plt.tight_layout()
    plt.savefig(output_path)
    print(f"Saved plot: {output_path}")

    # plot time by cluster type
    plt.figure(figsize=(10, 8.5))
    time_bins = np.linspace(0, 30, 15)
    for ctype in range(4):
        sel = clusterType_to_plot == ctype
        data = ak.flatten(time_to_plot[sel])
        counts, bin_edges = np.histogram(data, bins=time_bins)
        bin_centers = (bin_edges[:-1] + bin_edges[1:]) / 2
        errors = np.sqrt(counts)
        # Plot histogram and get the color
        hist_plot = plt.hist(data, bins=time_bins, histtype='step', label=cluster_labels[ctype], linewidth=1.5, alpha=0.8)
        hist_color = hist_plot[2][0].get_edgecolor()
        # Plot error bars with matching color and central marker
        plt.errorbar(bin_centers, counts, yerr=errors, fmt='o', markersize=3, elinewidth=1.2, capsize=2.5, capthick=1.2, color=hist_color, alpha=0.9)
    plt.xlabel('MergedCluster Time [ns]')
    plt.ylabel('Entries')
    plt.legend()
    # add cms label
    hep.cms.label("Preliminary", data=False, com = 13.6)
    # save figure
    output_path = output_dir / "time_byClusterType_reallyMergedOnly.png"
    plt.tight_layout()
    plt.savefig(output_path)
    print(f"Saved plot: {output_path}")

    # plot nClusters by cluster type
    plt.figure(figsize=(10, 8.5))
    ncluster_to_plot = nClusters[mask][good_evts]
    ncluster_bins = np.arange(-0.5, 10.5, 1)
    for ctype in range(4):
        sel = clusterType_to_plot == ctype
        data = ak.flatten(ncluster_to_plot[sel])
        counts, bin_edges = np.histogram(data, bins=ncluster_bins)
        bin_centers = (bin_edges[:-1] + bin_edges[1:]) / 2
        errors = np.sqrt(counts)
        # Plot histogram and get the color
        hist_plot = plt.hist(data, bins=ncluster_bins, histtype='step', label=cluster_labels[ctype], linewidth=1.5, alpha=0.8)
        hist_color = hist_plot[2][0].get_edgecolor()
        # Plot error bars with matching color and central marker
        plt.errorbar(bin_centers, counts, yerr=errors, fmt='o', markersize=3, elinewidth=1.2, capsize=2.5, capthick=1.2, color=hist_color, alpha=0.9)
    plt.xlabel('Number of clusters in Merged Clusters [ns]')
    plt.ylabel('Entries')
    plt.legend()
    # set y log scale
    plt.yscale('log')
    # add cms label
    hep.cms.label("Preliminary", data=False, com = 13.6)
    # save figure
    output_path = output_dir / "nClusters_byClusterType_reallyMergedOnly.png"
    plt.tight_layout()
    plt.savefig(output_path)
    print(f"Saved plot: {output_path}")

    # -----------------------------------------
    # -------- CLUSTER SIZE BY TYPE -----------
    # -----------------------------------------

    # plot cluster size in #crystals per cluster type
    mask = nClusters > 1
    good_evts = ak.sum(mask, axis = 1) > 0

    hitCols = arrays["simmc_hitCols_perCluster"]
    ieta = arrays["simmc_ieta_perCluster"]

    min_col = ak.min(ak.flatten(16*ieta + hitCols, axis = 3), axis = 2)
    max_col = ak.max(ak.flatten(16*ieta + hitCols, axis = 3), axis = 2)
    cluster_len = max_col - min_col + 1
    # apply mask and only consider merged clusters with >=1 clusters
    cluster_len = cluster_len[mask][good_evts]

    # finally, plot cluster size by cluster type
    plt.figure(figsize=(10, 8.5))
    size_bins = np.arange(0.5, 13.5, 1)
    for ctype in range(4):
        sel = clusterType_to_plot == ctype
        data = ak.flatten(cluster_len[sel])
        counts, bin_edges = np.histogram(data, bins=size_bins)
        bin_centers = (bin_edges[:-1] + bin_edges[1:]) / 2
        errors = np.sqrt(counts)
        # Plot histogram and get the color
        hist_plot = plt.hist(data, bins=size_bins, histtype='step', label=cluster_labels[ctype], linewidth=1.5, alpha=0.8)
        hist_color = hist_plot[2][0].get_edgecolor()
        # Plot error bars with matching color and central marker
        plt.errorbar(bin_centers, counts, yerr=errors, fmt='o', markersize=3, elinewidth=1.2, capsize=2.5, capthick=1.2, color=hist_color, alpha=0.9)

    plt.xlabel('MergedCluster Size [#crystal range]')
    plt.ylabel('Entries')
    plt.legend()
    plt.yscale('log')
    # add cms label
    hep.cms.label("Preliminary", data=False, com = 13.6)
    # save figure
    output_path = output_dir / "sizeInCrystals_byClusterType_reallyMergedOnly.png"
    plt.tight_layout()
    plt.savefig(output_path)
    print(f"Saved plot: {output_path}")

    # ---------------------------------

    # Also split by NUMBER OF MODULES
    # (=> count unique numbers in ieta)

    ### 1 module
    plt.figure(figsize=(10, 8.5))
    size_bins = np.arange(0.5, 13.5, 1)
    nmodules = arrays["simmc_nModules"][mask][good_evts]

    for ctype in range(4):
        sel = clusterType_to_plot == ctype
        sel2 = nmodules == 1
        data = ak.flatten(cluster_len[sel * sel2])
        counts, bin_edges = np.histogram(data, bins=size_bins)
        bin_centers = (bin_edges[:-1] + bin_edges[1:]) / 2
        errors = np.sqrt(counts)
        # Plot histogram and get the color
        hist_plot = plt.hist(data, bins=size_bins, histtype='step', label=cluster_labels[ctype], linewidth=1.5, alpha=0.8)
        hist_color = hist_plot[2][0].get_edgecolor()
        # Plot error bars with matching color and central marker
        plt.errorbar(bin_centers, counts, yerr=errors, fmt='o', markersize=3, elinewidth=1.2, capsize=2.5, capthick=1.2, color=hist_color, alpha=0.9)

    plt.xlabel('MergedCluster Size [#crystal range]')
    plt.ylabel('Entries')
    plt.legend()
    plt.yscale('log')
    # add cms label
    hep.cms.label("Preliminary", data=False, com = 13.6)
    # save figure
    output_path = output_dir / "sizeInCrystals_byClusterType_reallyMergedOnly_1module.png"
    plt.tight_layout()
    plt.savefig(output_path)
    print(f"Saved plot: {output_path}")

    ### >=2 modules
    plt.figure(figsize=(10, 8.5))
    size_bins = np.arange(0.5, 13.5, 1)
    nmodules = arrays["simmc_nModules"][mask][good_evts]

    for ctype in range(4):
        sel = clusterType_to_plot == ctype
        sel2 = nmodules >= 2
        data = ak.flatten(cluster_len[sel * sel2])
        counts, bin_edges = np.histogram(data, bins=size_bins)
        bin_centers = (bin_edges[:-1] + bin_edges[1:]) / 2
        errors = np.sqrt(counts)
        # Plot histogram and get the color
        hist_plot = plt.hist(data, bins=size_bins, histtype='step', label=cluster_labels[ctype], linewidth=1.5, alpha=0.8)
        hist_color = hist_plot[2][0].get_edgecolor()
        # Plot error bars with matching color and central marker
        plt.errorbar(bin_centers, counts, yerr=errors, fmt='o', markersize=3, elinewidth=1.2, capsize=2.5, capthick=1.2, color=hist_color, alpha=0.9)

    plt.xlabel('MergedCluster Size [#crystal range]')
    plt.ylabel('Entries')
    plt.legend()
    plt.yscale('log')
    # add cms label
    hep.cms.label("Preliminary", data=False, com = 13.6)
    # save figure
    output_path = output_dir / "sizeInCrystals_byClusterType_reallyMergedOnly_2+modules.png"
    plt.tight_layout()
    plt.savefig(output_path)
    print(f"Saved plot: {output_path}")    

    # ---------------------------------

    # For comparison, plot cluster size for 1-cluster MCs
    mask = nClusters == 1
    good_evts = ak.sum(mask, axis = 1) > 0
    cluster_len = max_col - min_col + 1
    cluster_len = cluster_len[mask][good_evts]
    clusterType_to_plot = ak.firsts(clusterTypes[mask][good_evts], axis = 2)

    plt.figure(figsize=(10, 8.5))
    size_bins = np.arange(0.5, 13.5, 1)
    for ctype in range(4):
        sel = clusterType_to_plot == ctype
        data = ak.flatten(cluster_len[sel])
        counts, bin_edges = np.histogram(data, bins=size_bins)
        bin_centers = (bin_edges[:-1] + bin_edges[1:]) / 2
        errors = np.sqrt(counts)
        # Plot histogram and get the color
        hist_plot = plt.hist(data, bins=size_bins, histtype='step', label=cluster_labels[ctype], linewidth=1.5, alpha=0.8)
        hist_color = hist_plot[2][0].get_edgecolor()
        # Plot error bars with matching color and central marker
        plt.errorbar(bin_centers, counts, yerr=errors, fmt='o', markersize=3, elinewidth=1.2, capsize=2.5, capthick=1.2, color=hist_color, alpha=0.9)

    plt.xlabel('MergedCluster Size [#crystal range]')
    plt.ylabel('Entries')
    plt.legend()
    plt.yscale('log')
    # add cms label
    hep.cms.label("Preliminary", data=False, com = 13.6)
    # save figure
    output_path = output_dir / "sizeInCrystals_byClusterType_1cluster.png"
    plt.tight_layout()
    plt.savefig(output_path)
    print(f"Saved plot: {output_path}")

    print(f"\nAll plots saved to {output_dir}")
    return 0

if __name__ == "__main__":
    exit(main())
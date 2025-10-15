import numpy as np

PLOT_DEFINITIONS = [
    # 1D histograms
    {
        'branches': ['simmc_n'],
        'name': 'simmc_multiplicity',
        'xlabel': 'Number of MergedClusters',
        'ylabel': 'Events',
        'binning': np.arange(-0.5, 50.5, 1),
    },
    {
        'branches': ['simmc_energy'],
        'name': 'simmc_energy',
        'xlabel': 'MergedCluster Energy [MeV]',
        'ylabel': 'MergedClusters',
        'binning': np.linspace(0, 5, 50),
        'cuts': {'simmc_nClusters': lambda x: x > 0}  # Only mergedclusters with >0 clusters
    },
    {
        'branches': ['simmc_energy'],
        'name': 'simmc_logEnergy',
        'xlabel': 'MergedCluster Energy [MeV]',
        'ylabel': 'MergedClusters',
        'binning': np.logspace(-2, 3, 50),
        'logx_scale' : True,
        'cuts': {'simmc_nClusters': lambda x: x > 0}  # Only mergedclusters with >0 clusters
    },
    {
        'branches': ['simmc_energy'],
        'name': 'simmc_energy_log',
        'xlabel': 'MergedCluster Energy [MeV]',
        'ylabel': 'MergedClusters',
        'binning': np.linspace(0, 100, 50),
        'logy_scale': True,
        'cuts': {'simmc_nClusters': lambda x: x > 0}  # Only mergedclusters with >0 clusters
    },
    {
        'branches': ['simmc_time'],
        'name': 'simmc_time', 
        'xlabel': 'MergedCluster Time [ns]',
        'ylabel': 'MergedClusters',
        'binning': np.linspace(0, 30, 51),
        'cuts': {'simmc_nClusters': lambda x: x > 0}  # Only mergedclusters with >0 clusters
    },
    {
        'branches': ['simmc_x'],
        'name': 'simmc_x',
        'xlabel': 'MergedCluster X [mm]',
        'ylabel': 'MergedClusters', 
        'binning': np.linspace(-120, 120, 51),
        'cuts': {'simmc_nClusters': lambda x: x > 0}  # Only mergedclusters with >0 clusters
    },
    {
        'branches': ['simmc_y'],
        'name': 'simmc_y',
        'xlabel': 'MergedCluster Y [mm]',
        'ylabel': 'MergedClusters',
        'binning': np.linspace(-120, 120, 51),
        'cuts': {'simmc_nClusters': lambda x: x > 0}  # Only mergedclusters with >0 clusters
    },
    {
        'branches': ['simmc_nClusters'],
        'name': 'simmc_nClusters',
        'xlabel': 'Number of Clusters in MergedCluster',
        'ylabel': 'MergedClusters',
        'binning': np.arange(-0.5, 100.5, 1),
        'logy_scale': True
    },

    # 1D histograms -- per-cluster quantities
    {
        'branches': ['simmc_energy_perCluster'],
        'name': 'simmc_energy_perCluster',
        'xlabel': 'Cluster energy [MeV]',
        'ylabel': 'Clusters',
        'binning': np.logspace(-3, 3, 100),
        'type_breakdown': True,
        'logx_scale': True,
        'cuts': {'simmc_nClusters': lambda x: x > 0}  # Only clusters in mergedclusters with >0 clusters
    },

    {
        'branches': ['simmc_energy_perCluster'],
        'name': 'simmc_energy_perCluster_log',
        'xlabel': 'Cluster energy [MeV]',
        'ylabel': 'Clusters',
        'binning': np.logspace(-3, 3, 100),
        'type_breakdown': True,
        'logx_scale': True,
        'logy_scale': True,
        'cuts': {'simmc_nClusters': lambda x: x > 0}  # Only clusters in mergedclusters with >0 clusters
    },

    {
        'branches': ['simmc_time_perCluster'],
        'name': 'simmc_time_perCluster',
        'xlabel': 'Cluster time [ns]',
        'ylabel': 'Clusters',
        'binning': np.linspace(0, 30, 51),
        'type_breakdown': True,
        'cuts': {'simmc_nClusters': lambda x: x > 0}  # Only clusters in mergedclusters with >0 clusters
    },

    {
        'branches': ['simmc_time_perCluster'],
        'name': 'simmc_time_perCluster_log',
        'xlabel': 'Cluster time [ns]',
        'ylabel': 'Clusters',
        'binning': np.linspace(0, 30, 51),
        'type_breakdown': True,
        'logy_scale': True,
        'cuts': {'simmc_nClusters': lambda x: x > 0}  # Only clusters in mergedclusters with >0 clusters
    },

    {
        'branches' : ['simmc_clusterType'],
        'name' : 'simmc_clusterType',
        'xlabel' : 'Cluster type',
        'ylabel' : 'Clusters',
        'binning' : np.arange(-0.5, 4.5, 1),
        'logy_scale': True,
        'cuts': {'simmc_nClusters': lambda x: x > 0}  #
    },

    # 2D histograms
    {
        'branches': ['simmc_x', 'simmc_y'],
        'name': 'simmc_position_xy',
        'xlabel': 'MergedCluster X [mm]',
        'ylabel': 'MergedCluster Y [mm]',
        'binning': [np.linspace(-120, 120, 101), np.linspace(-120, 120, 101)],
        'cuts': {'simmc_nClusters': lambda x: x > 0}  # Only mergedclusters with >0 clusters
    },
    {
        'branches': ['simmc_time', 'simmc_energy'],
        'name': 'simmc_energy_vs_time',
        'xlabel': 'MergedCluster Time [ns]',
        'ylabel': 'MergedCluster Energy [MeV]',
        'binning': [np.linspace(4, 27, 31), np.linspace(0, 1, 21)],
        'cuts': {'simmc_nClusters': lambda x: x > 0}  # Only mergedclusters with >0 clusters
    },

    {
        'branches': ['simmc_nClusters', 'simmc_energy'],
        'name': 'simmc_energy_vs_nClusters',
        'xlabel': 'Number of Clusters',
        'ylabel': 'MergedCluster Energy [MeV]',
        'binning': [np.arange(-0.5, 10.5, 1), np.linspace(0, 3, 21)],
        'logy_scale': True,
        'cuts': {'simmc_nClusters': lambda x: x > 0}  # Only mergedclusters with >0 clusters
    },

    {
        'branches': ['simmc_nClusters', 'simmc_primary_et'],
        'name': 'simmc_primaryPt_vs_nClusters',
        'xlabel': 'Number of Clusters',
        'ylabel': 'Primary particle pT [GeV]',
        'binning': [np.arange(-0.5, 20.5, 1), np.linspace(0, 11, 21)],
        'cuts': {'simmc_nClusters': lambda x: x > 0}  # Only mergedclusters with >0 clusters
    },

    {
        'branches': ['simmc_nClusters', 'simmc_primary_energy'],
        'name': 'simmc_primaryEnergy_vs_nClusters',
        'xlabel': 'Number of Clusters',
        'ylabel': 'Primary particle energy [GeV]',
        'binning': [np.arange(-0.5, 20.5, 1), np.linspace(0, 40, 21)],
        'cuts': {'simmc_nClusters': lambda x: x > 0}  # Only mergedclusters with >0 clusters
    },

    {
        'branches': ['simmc_energy', 'simmc_primary_et'],
        'name': 'simmc_primaryPt_vs_energy',
        'xlabel': 'Supercluster energy [MeV]',
        'ylabel': 'Primary particle pT [GeV]',
        'binning': [np.linspace(0, 3, 21), np.linspace(0, 11, 21)],
        'cuts': {'simmc_nClusters': lambda x: x > 0}  # Only mergedclusters with >0 clusters
    },

    {
        'branches': ['simmc_energy', 'simmc_primary_energy'],
        'name': 'simmc_primaryEnergy_vs_energy',
        'xlabel': 'SC energy [MeV]',
        'ylabel': 'Primary particle energy [GeV]',
        'binning': [np.linspace(0, 5, 50), np.linspace(0, 40, 21)],
        'cuts': {'simmc_nClusters': lambda x: x > 0}  # Only mergedclusters with >0 clusters
    },

]
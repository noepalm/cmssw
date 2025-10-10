import numpy as np

PLOT_DEFINITIONS = [
    # 1D histograms
    {
        'branches': ['sc_n'],
        'name': 'sc_multiplicity',
        'xlabel': 'Number of SuperClusters',
        'ylabel': 'Events',
        'binning': np.arange(-0.5, 50.5, 1),
    },
    {
        'branches': ['sc_energy'],
        'name': 'sc_energy',
        'xlabel': 'SuperCluster Energy [MeV]',
        'ylabel': 'SuperClusters',
        'binning': np.linspace(0, 5, 50),
        'cuts': {'sc_nClusters': lambda x: x > 0}  # Only superclusters with >0 clusters
    },
    {
        'branches': ['sc_energy'],
        'name': 'sc_logEnergy',
        'xlabel': 'SuperCluster Energy [MeV]',
        'ylabel': 'SuperClusters',
        'binning': np.logspace(-2, 3, 50),
        'logx_scale' : True,
        'cuts': {'sc_nClusters': lambda x: x > 0}  # Only superclusters with >0 clusters
    },
    {
        'branches': ['sc_energy'],
        'name': 'sc_energy_log',
        'xlabel': 'SuperCluster Energy [MeV]',
        'ylabel': 'SuperClusters',
        'binning': np.linspace(0, 100, 50),
        'logy_scale': True,
        'cuts': {'sc_nClusters': lambda x: x > 0}  # Only superclusters with >0 clusters
    },
    {
        'branches': ['sc_time'],
        'name': 'sc_time', 
        'xlabel': 'SuperCluster Time [ns]',
        'ylabel': 'SuperClusters',
        'binning': np.linspace(0, 30, 51),
        'cuts': {'sc_nClusters': lambda x: x > 0}  # Only superclusters with >0 clusters
    },
    {
        'branches': ['sc_x'],
        'name': 'sc_x',
        'xlabel': 'SuperCluster X [mm]',
        'ylabel': 'SuperClusters', 
        'binning': np.linspace(-120, 120, 51),
        'cuts': {'sc_nClusters': lambda x: x > 0}  # Only superclusters with >0 clusters
    },
    {
        'branches': ['sc_y'],
        'name': 'sc_y',
        'xlabel': 'SuperCluster Y [mm]',
        'ylabel': 'SuperClusters',
        'binning': np.linspace(-120, 120, 51),
        'cuts': {'sc_nClusters': lambda x: x > 0}  # Only superclusters with >0 clusters
    },
    {
        'branches': ['sc_nClusters'],
        'name': 'sc_nClusters',
        'xlabel': 'Number of Clusters in SuperCluster',
        'ylabel': 'SuperClusters',
        'binning': np.arange(-0.5, 100.5, 1),
        'logy_scale': True
    },

    # 1D histograms -- per-cluster quantities
    {
        'branches': ['sc_energy_perCluster'],
        'name': 'sc_energy_perCluster',
        'xlabel': 'Cluster energy [MeV]',
        'ylabel': 'Clusters',
        'binning': np.logspace(-3, 3, 100),
        'type_breakdown': True,
        'logx_scale': True,
        'cuts': {'sc_nClusters': lambda x: x > 0}  # Only clusters in superclusters with >0 clusters
    },

    {
        'branches': ['sc_energy_perCluster'],
        'name': 'sc_energy_perCluster_log',
        'xlabel': 'Cluster energy [MeV]',
        'ylabel': 'Clusters',
        'binning': np.logspace(-3, 3, 100),
        'type_breakdown': True,
        'logx_scale': True,
        'logy_scale': True,
        'cuts': {'sc_nClusters': lambda x: x > 0}  # Only clusters in superclusters with >0 clusters
    },

    {
        'branches': ['sc_time_perCluster'],
        'name': 'sc_time_perCluster',
        'xlabel': 'Cluster time [ns]',
        'ylabel': 'Clusters',
        'binning': np.linspace(0, 30, 51),
        'type_breakdown': True,
        'cuts': {'sc_nClusters': lambda x: x > 0}  # Only clusters in superclusters with >0 clusters
    },

    {
        'branches': ['sc_time_perCluster'],
        'name': 'sc_time_perCluster_log',
        'xlabel': 'Cluster time [ns]',
        'ylabel': 'Clusters',
        'binning': np.linspace(0, 30, 51),
        'type_breakdown': True,
        'logy_scale': True,
        'cuts': {'sc_nClusters': lambda x: x > 0}  # Only clusters in superclusters with >0 clusters
    },

    {
        'branches' : ['sc_clusterType'],
        'name' : 'sc_clusterType',
        'xlabel' : 'Cluster type',
        'ylabel' : 'Clusters',
        'binning' : np.arange(-0.5, 4.5, 1),
        'logy_scale': True,
        'cuts': {'sc_nClusters': lambda x: x > 0}  #
    },

    # 2D histograms
    {
        'branches': ['sc_x', 'sc_y'],
        'name': 'sc_position_xy',
        'xlabel': 'SuperCluster X [mm]',
        'ylabel': 'SuperCluster Y [mm]',
        'binning': [np.linspace(-120, 120, 101), np.linspace(-120, 120, 101)],
        'cuts': {'sc_nClusters': lambda x: x > 0}  # Only superclusters with >0 clusters
    },
    {
        'branches': ['sc_time', 'sc_energy'],
        'name': 'sc_energy_vs_time',
        'xlabel': 'SuperCluster Time [ns]',
        'ylabel': 'SuperCluster Energy [MeV]',
        'binning': [np.linspace(4, 27, 31), np.linspace(0, 1, 21)],
        'cuts': {'sc_nClusters': lambda x: x > 0}  # Only superclusters with >0 clusters
    },

    {
        'branches': ['sc_nClusters', 'sc_energy'],
        'name': 'sc_energy_vs_nClusters',
        'xlabel': 'Number of Clusters',
        'ylabel': 'SuperCluster Energy [MeV]',
        'binning': [np.arange(-0.5, 10.5, 1), np.linspace(0, 3, 21)],
        'logy_scale': True,
        'cuts': {'sc_nClusters': lambda x: x > 0}  # Only superclusters with >0 clusters
    },

    {
        'branches': ['sc_nClusters', 'sc_primary_et'],
        'name': 'sc_primaryPt_vs_nClusters',
        'xlabel': 'Number of Clusters',
        'ylabel': 'Primary particle pT [GeV]',
        'binning': [np.arange(-0.5, 20.5, 1), np.linspace(0, 11, 21)],
        'cuts': {'sc_nClusters': lambda x: x > 0}  # Only superclusters with >0 clusters
    },

    {
        'branches': ['sc_nClusters', 'sc_primary_energy'],
        'name': 'sc_primaryEnergy_vs_nClusters',
        'xlabel': 'Number of Clusters',
        'ylabel': 'Primary particle energy [GeV]',
        'binning': [np.arange(-0.5, 20.5, 1), np.linspace(0, 40, 21)],
        'cuts': {'sc_nClusters': lambda x: x > 0}  # Only superclusters with >0 clusters
    },

    {
        'branches': ['sc_energy', 'sc_primary_et'],
        'name': 'sc_primaryPt_vs_energy',
        'xlabel': 'Supercluster energy [MeV]',
        'ylabel': 'Primary particle pT [GeV]',
        'binning': [np.linspace(0, 3, 21), np.linspace(0, 11, 21)],
        'cuts': {'sc_nClusters': lambda x: x > 0}  # Only superclusters with >0 clusters
    },

    {
        'branches': ['sc_energy', 'sc_primary_energy'],
        'name': 'sc_primaryEnergy_vs_energy',
        'xlabel': 'SC energy [MeV]',
        'ylabel': 'Primary particle energy [GeV]',
        'binning': [np.linspace(0, 5, 50), np.linspace(0, 40, 21)],
        'cuts': {'sc_nClusters': lambda x: x > 0}  # Only superclusters with >0 clusters
    },

]
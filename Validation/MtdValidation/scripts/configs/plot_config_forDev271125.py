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
        'branches': ['simmc_eta'],
        'name': 'simmc_eta',
        'xlabel': r'$\eta$ of MergedClusters',
        'ylabel': 'Events',
        'binning': np.linspace(-1.5, 1.5, 20),
        'cuts': {'simmc_nClusters': lambda x: x > 0}
    },
    {
        'branches': ['simmc_eta'],
        'name': 'simmc_eta_reallyMergedOnly',
        'xlabel': r'$\eta$ of MergedClusters',
        'ylabel': 'Events',
        'binning': np.linspace(-1.5, 1.5, 20),
        'cuts': {'simmc_nClusters': lambda x: x > 1}
    },
    {
        'branches': ['simmc_energy'],
        'name': 'simmc_logEnergy_reallyMergedOnly',
        'xlabel': 'MergedCluster Energy [MeV]',
        'ylabel': 'MergedClusters',
        'binning': np.logspace(-3, 3, 50),
        'logx_scale' : True,
        'cuts': {'simmc_nClusters': lambda x: x > 1}
    },
    # {
    #     'branches': ['simmc_energy'],
    #     'name': 'simmc_logEnergy_reallyMergedOnly_perType',
    #     'xlabel': 'MergedCluster Energy [MeV]',
    #     'ylabel': 'MergedClusters',
    #     'binning': np.logspace(-3, 3, 15),
    #     'logx_scale' : True,
    #     'type_breakdown': True,
    #     'cuts': {'simmc_nClusters': lambda x: x > 1}
    # },    
    {
        'branches': ['simmc_energy_perCluster'],
        'name': 'simmc_logEnergy_perCluster',
        'xlabel': 'MergedCluster Energy [MeV]',
        'ylabel': 'MergedClusters',
        'binning': np.logspace(-3, 3, 100),
        'logx_scale' : True,
        'type_breakdown': True,
        'cuts': {'simmc_nClusters': lambda x: x > 0}
    },
    {
        'branches': ['simmc_energy_perCluster'],
        'name': 'simmc_logEnergy_perCluster',
        'xlabel': 'MergedCluster Energy [MeV]',
        'ylabel': 'MergedClusters',
        'binning': np.logspace(-3, 3, 100),
        'logx_scale' : True,
        'type_breakdown': True,
        'cuts': {'simmc_nClusters': lambda x: x > 0}
    },
    {
        'branches' : ['simmc_clusterType'],
        'name' : 'simmc_clusterType',
        'xlabel' : 'Cluster type',
        'ylabel' : 'Clusters',
        'binning' : np.arange(-0.5, 4.5, 1),
        'logy_scale': True,
        'plot_numbers': True,
        'cuts': {'simmc_nClusters': lambda x: x > 0}
    },
    {
        'branches': ['simmc_time'],
        'name': 'simmc_time_reallyMergedOnly', 
        'xlabel': 'MergedCluster Time [ns]',
        'ylabel': 'MergedClusters',
        'binning': np.linspace(0, 30, 51),
        'cuts': {'simmc_nClusters': lambda x: x > 1}
    },
    # {
    #     'branches': ['simmc_time'],
    #     'name': 'simmc_time_reallyMergedOnly_perType', 
    #     'xlabel': 'MergedCluster Time [ns]',
    #     'ylabel': 'MergedClusters',
    #     'binning': np.linspace(0, 30, 15),
    #     'type_breakdown': True,
    #     'cuts': {'simmc_nClusters': lambda x: x > 1}
    # },    
    {
        'branches': ['simmc_time_perCluster'],
        'name': 'simmc_time_perCluster', 
        'xlabel': 'MergedCluster Time [ns]',
        'ylabel': 'MergedClusters',
        'type_breakdown': True,
        'binning': np.linspace(0, 30, 51),
        'cuts': {'simmc_nClusters': lambda x: x > 0}
    },
    {
        'branches': ['simmc_time_perCluster'],
        'name': 'simmc_time_perCluster_logy', 
        'xlabel': 'MergedCluster Time [ns]',
        'ylabel': 'MergedClusters',
        'type_breakdown': True,
        'binning': np.linspace(0, 30, 51),
        'logy_scale': True,
        'cuts': {'simmc_nClusters': lambda x: x > 0}
    },    
    {
        'branches': ['simmc_nClusters'],
        'name': 'simmc_nClusters',
        'xlabel': 'Number of Clusters in MergedCluster',
        'ylabel': 'MergedClusters',
        'binning': np.arange(-0.5, 10.5, 1),
        'logy_scale': True,
        'plot_numbers': True,
        'cuts': {'simmc_nClusters': lambda x: x > 0}
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
        'binning': [np.linspace(4, 27, 31), np.linspace(-2, np.log(20), 21)],
        'logy_scale': True,
        'cuts': {'simmc_nClusters': lambda x: x > 0}  # Only mergedclusters with >0 clusters
    },
]

PLOT_DEFINITIONS.append({
    'branches': ['simmc_time_perCluster'],
    'name': 'simmc_time_perCluster',
    'xlabel': 'Cluster time [ns]',
    'ylabel': 'Entries',
    'binning': np.linspace(0, 30, 51),
    'cuts': {'simmc_nClusters': lambda x: x > 0},
    'overlap': [
        {
            'branches': ['simmc_time'],
            'name': 'simmc_time',
            'binning': np.linspace(0, 30, 51),
            'cuts': {'simmc_nClusters': lambda x: x > 0}
        }
    ],
})

PLOT_DEFINITIONS.append({
    'branches': ['simmc_energy_perCluster'],
    'name': 'simmc_energy_sanityCheck',
    'xlabel': 'Cluster energy [ns]',
    'ylabel': 'Entries',
    'binning': np.logspace(-3, 3, 100),
    'logx_scale' : True,
    'cuts': {'simmc_nClusters': lambda x: x > 0},
    'overlap': [
        {
            'branches': ['simmc_energy'],
            'name': 'simmc_energy',
            'binning': np.logspace(-3, 3, 100),
            'cuts': {'simmc_nClusters': lambda x: x > 0}
        }
    ],
})
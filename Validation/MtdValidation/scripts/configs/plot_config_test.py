import numpy as np

PLOT_DEFINITIONS = [
    {
        'branches': ['simmc_time_perCluster'],
        'name': 'simmc_time_perCluster_logy', 
        'xlabel': 'MergedCluster Time [ns]',
        'ylabel': 'MergedClusters',
        'type_breakdown': True,
        'binning': np.linspace(0, 30, 51),
        'logy_scale': True,
        'cuts': {'simmc_nClusters': lambda x: x > 0},
    },
    {
        'branches': ['simmc_time_perCluster'],
        'name': 'simmc_time_perCluster_logy_eta<0.1', 
        'xlabel': 'MergedCluster Time [ns]',
        'ylabel': 'MergedClusters',
        'type_breakdown': True,
        'binning': np.linspace(0, 30, 51),
        'logy_scale': True,
        'cuts': {'simmc_nClusters': lambda x: x > 0, 'simmc_eta' : lambda x: abs(x) < 0.1},
    },
    {
        'branches': ['simmc_time_perCluster'],
        'name': 'simmc_time_perCluster_logy_0.1<eta<0.2', 
        'xlabel': 'MergedCluster Time [ns]',
        'ylabel': 'MergedClusters',
        'type_breakdown': True,
        'binning': np.linspace(0, 30, 51),
        'logy_scale': True,
        'cuts': {'simmc_nClusters': lambda x: x > 0, 'simmc_eta' : lambda x: (abs(x) < 0.2) * (abs(x) >= 0.1)},
    },
    {
        'branches': ['simmc_time_perCluster'],
        'name': 'simmc_time_perCluster_logy_0.2<eta<0.3', 
        'xlabel': 'MergedCluster Time [ns]',
        'ylabel': 'MergedClusters',
        'type_breakdown': True,
        'binning': np.linspace(0, 30, 51),
        'logy_scale': True,
        'cuts': {'simmc_nClusters': lambda x: x > 0, 'simmc_eta' : lambda x: (abs(x) < 0.3) * (abs(x) >= 0.2)},
    },
    {
        'branches': ['simmc_time_perCluster'],
        'name': 'simmc_time_perCluster_logy_eta>1.3', 
        'xlabel': 'MergedCluster Time [ns]',
        'ylabel': 'MergedClusters',
        'type_breakdown': True,
        'binning': np.linspace(0, 30, 51),
        'logy_scale': True,
        'cuts': {'simmc_nClusters': lambda x: x > 0, 'simmc_eta' : lambda x: (abs(x) > 1.3) * (abs(x) < 1.5)},
    },


]
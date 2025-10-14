import numpy as np

from configs.plot_config_history import PLOT_DEFINITIONS as PLOT_DEFINITIONS_HISTORY

new_binnings = {
    "simmc_energy_log" : np.linspace(0, 50, 50),
    "simmc_logEnergy" : np.logspace(-3, 3, 100),
    'simmc_nClusters' : np.arange(-0.5, 10.5, 1),
    'simmc_energy_vs_nClusters' : [np.arange(-0.5, 3.5, 1), np.linspace(0, 1, 20)],
    'simmc_primaryPt_vs_nClusters' : [np.arange(-0.5, 3.5, 1), np.linspace(0, 11, 20)],
    'simmc_primaryEnergy_vs_nClusters' : [np.arange(-0.5, 3.5, 1), np.linspace(0, 40, 20)],
    'simmc_primaryEnergy_vs_energy' : [np.linspace(0, 1, 20), np.linspace(0, 40, 20)],
    'simmc_primaryPt_vs_energy' : [np.linspace(0, 1, 20), np.linspace(0, 11, 21)],
}

for name, binning in new_binnings.items():
    # find plot entry with corresponding name and update its binning
    # NOTE: will throw StopIteration error if name not found
    next(filter(lambda x: x.get("name", None) == name, PLOT_DEFINITIONS_HISTORY))['binning'] = binning

PLOT_DEFINITIONS = PLOT_DEFINITIONS_HISTORY

# [
#     # 1D histograms
#     {
#         'branches': ['simmc_n'],
#         'name': 'simmc_multiplicity',
#         'xlabel': 'Number of MergedClusters',
#         'ylabel': 'Events',
#         'binning': np.arange(-0.5, 50.5, 1),
#         'log_scale': False
#     },
#     {
#         'branches': ['simmc_energy'],
#         'name': 'simmc_energy',
#         'xlabel': 'MergedCluster Energy [MeV]',
#         'ylabel': 'MergedClusters',
#         'binning': np.linspace(0, 5, 50),
#         'log_scale': False,
#         'cuts': {'simmc_nClusters': lambda x: x > 0}  # Only mergedclusters with >0 clusters
#     },
#     {
#         'branches': ['simmc_energy'],
#         'name': 'simmc_energy',
#         'xlabel': 'MergedCluster Energy [MeV]',
#         'ylabel': 'MergedClusters',
#         'binning': np.linspace(0, 50, 50),
#         'log_scale': True,
#         'cuts': {'simmc_nClusters': lambda x: x > 0}  # Only mergedclusters with >0 clusters
#     },
#     {
#         'branches': ['simmc_time'],
#         'name': 'simmc_time', 
#         'xlabel': 'MergedCluster Time [ns]',
#         'ylabel': 'MergedClusters',
#         'binning': np.linspace(0, 30, 51),
#         'log_scale': False,
#         'cuts': {'simmc_nClusters': lambda x: x > 0}  # Only mergedclusters with >0 clusters
#     },
#     {
#         'branches': ['simmc_x'],
#         'name': 'simmc_x',
#         'xlabel': 'MergedCluster X [mm]',
#         'ylabel': 'MergedClusters', 
#         'binning': np.linspace(-120, 120, 51),
#         'log_scale': False,
#         'cuts': {'simmc_nClusters': lambda x: x > 0}  # Only mergedclusters with >0 clusters
#     },
#     {
#         'branches': ['simmc_y'],
#         'name': 'simmc_y',
#         'xlabel': 'MergedCluster Y [mm]',
#         'ylabel': 'MergedClusters',
#         'binning': np.linspace(-120, 120, 51),
#         'log_scale': False,
#         'cuts': {'simmc_nClusters': lambda x: x > 0}  # Only mergedclusters with >0 clusters
#     },
#     {
#         'branches': ['simmc_nClusters'],
#         'name': 'simmc_nClusters',
#         'xlabel': 'Number of Clusters in MergedCluster',
#         'ylabel': 'MergedClusters',
#         'binning': np.arange(-0.5, 10.5, 1),
#         'log_scale': True
#     },

#     # 2D histograms
#     {
#         'branches': ['simmc_x', 'simmc_y'],
#         'name': 'simmc_position_xy',
#         'xlabel': 'MergedCluster X [mm]',
#         'ylabel': 'MergedCluster Y [mm]',
#         'binning': [np.linspace(-120, 120, 101), np.linspace(-120, 120, 101)],
#         'log_scale': False,
#         'cuts': {'simmc_nClusters': lambda x: x > 0}  # Only mergedclusters with >0 clusters
#     },
#     {
#         'branches': ['simmc_time', 'simmc_energy'],
#         'name': 'simmc_energy_vs_time',
#         'xlabel': 'MergedCluster Time [ns]',
#         'ylabel': 'MergedCluster Energy [MeV]',
#         'binning': [np.linspace(4, 27, 31), np.linspace(0, 1, 21)],
#         'log_scale': False,
#         'cuts': {'simmc_nClusters': lambda x: x > 0}  # Only mergedclusters with >0 clusters
#     },

#     {
#         'branches': ['simmc_nClusters', 'simmc_energy'],
#         'name': 'simmc_energy_vs_nClusters',
#         'xlabel': 'Number of Clusters',
#         'ylabel': 'MergedCluster Energy [MeV]',
#         'binning': [np.arange(-0.5, 3.5, 1), np.linspace(0, 1, 21)],
#         'log_scale': True,
#         'cuts': {'simmc_nClusters': lambda x: x > 0}  # Only mergedclusters with >0 clusters
#     }


#     {
#         'branches': ['simmc_nClusters', 'sd_primary_et'],
#         'name': 'simmc_primaryPt_vs_nClusters',
#         'xlabel': 'Number of Clusters',
#         'ylabel': 'Primary particle pT [GeV]',
#         'binning': [np.arange(-0.5, 3.5, 1), np.linspace(0, 11, 21)],
#         'cuts': {'simmc_nClusters': lambda x: x > 0}  # Only mergedclusters with >0 clusters
#     },

#     {
#         'branches': ['simmc_nClusters', 'sd_primary_energy'],
#         'name': 'simmc_primaryEnergy_vs_nClusters',
#         'xlabel': 'Number of Clusters',
#         'ylabel': 'Primary particle energy [GeV]',
#         'binning': [np.arange(-0.5, 3.5, 1), np.linspace(0, 40, 21)],
#         'cuts': {'simmc_nClusters': lambda x: x > 0}  # Only mergedclusters with >0 clusters
#     },

#     {
#         'branches': ['simmc_energy', 'sd_primary_energy'],
#         'name': 'simmc_energy_vs_nClusters',
#         'xlabel': 'SC energy [MeV]',
#         'ylabel': 'Primary particle energy [GeV]',
#         'binning': [np.linspace(0, 5, 50), np.linspace(0, 40, 21)],
#         'cuts': {'simmc_nClusters': lambda x: x > 0}  # Only mergedclusters with >0 clusters
#     },

# ]
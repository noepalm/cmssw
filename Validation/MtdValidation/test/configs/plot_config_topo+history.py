import numpy as np

from configs.plot_config_history import PLOT_DEFINITIONS as PLOT_DEFINITIONS_HISTORY

new_binnings = {
    "sc_energy_log" : np.linspace(0, 50, 50),
    "sc_logEnergy" : np.logspace(-3, 3, 100),
    'sc_nClusters' : np.arange(-0.5, 10.5, 1),
    'sc_energy_vs_nClusters' : [np.arange(-0.5, 3.5, 1), np.linspace(0, 1, 20)],
    'sc_primaryPt_vs_nClusters' : [np.arange(-0.5, 3.5, 1), np.linspace(0, 11, 20)],
    'sc_primaryEnergy_vs_nClusters' : [np.arange(-0.5, 3.5, 1), np.linspace(0, 40, 20)],
    'sc_primaryEnergy_vs_energy' : [np.linspace(0, 1, 20), np.linspace(0, 40, 20)],
    'sc_primaryPt_vs_energy' : [np.linspace(0, 1, 20), np.linspace(0, 11, 21)],
}

for name, binning in new_binnings.items():
    # find plot entry with corresponding name and update its binning
    # NOTE: will throw StopIteration error if name not found
    next(filter(lambda x: x.get("name", None) == name, PLOT_DEFINITIONS_HISTORY))['binning'] = binning

PLOT_DEFINITIONS = PLOT_DEFINITIONS_HISTORY

# [
#     # 1D histograms
#     {
#         'branches': ['sc_n'],
#         'name': 'sc_multiplicity',
#         'xlabel': 'Number of SuperClusters',
#         'ylabel': 'Events',
#         'binning': np.arange(-0.5, 50.5, 1),
#         'log_scale': False
#     },
#     {
#         'branches': ['sc_energy'],
#         'name': 'sc_energy',
#         'xlabel': 'SuperCluster Energy [MeV]',
#         'ylabel': 'SuperClusters',
#         'binning': np.linspace(0, 5, 50),
#         'log_scale': False,
#         'cuts': {'sc_nClusters': lambda x: x > 0}  # Only superclusters with >0 clusters
#     },
#     {
#         'branches': ['sc_energy'],
#         'name': 'sc_energy',
#         'xlabel': 'SuperCluster Energy [MeV]',
#         'ylabel': 'SuperClusters',
#         'binning': np.linspace(0, 50, 50),
#         'log_scale': True,
#         'cuts': {'sc_nClusters': lambda x: x > 0}  # Only superclusters with >0 clusters
#     },
#     {
#         'branches': ['sc_time'],
#         'name': 'sc_time', 
#         'xlabel': 'SuperCluster Time [ns]',
#         'ylabel': 'SuperClusters',
#         'binning': np.linspace(0, 30, 51),
#         'log_scale': False,
#         'cuts': {'sc_nClusters': lambda x: x > 0}  # Only superclusters with >0 clusters
#     },
#     {
#         'branches': ['sc_x'],
#         'name': 'sc_x',
#         'xlabel': 'SuperCluster X [mm]',
#         'ylabel': 'SuperClusters', 
#         'binning': np.linspace(-120, 120, 51),
#         'log_scale': False,
#         'cuts': {'sc_nClusters': lambda x: x > 0}  # Only superclusters with >0 clusters
#     },
#     {
#         'branches': ['sc_y'],
#         'name': 'sc_y',
#         'xlabel': 'SuperCluster Y [mm]',
#         'ylabel': 'SuperClusters',
#         'binning': np.linspace(-120, 120, 51),
#         'log_scale': False,
#         'cuts': {'sc_nClusters': lambda x: x > 0}  # Only superclusters with >0 clusters
#     },
#     {
#         'branches': ['sc_nClusters'],
#         'name': 'sc_nClusters',
#         'xlabel': 'Number of Clusters in SuperCluster',
#         'ylabel': 'SuperClusters',
#         'binning': np.arange(-0.5, 10.5, 1),
#         'log_scale': True
#     },

#     # 2D histograms
#     {
#         'branches': ['sc_x', 'sc_y'],
#         'name': 'sc_position_xy',
#         'xlabel': 'SuperCluster X [mm]',
#         'ylabel': 'SuperCluster Y [mm]',
#         'binning': [np.linspace(-120, 120, 101), np.linspace(-120, 120, 101)],
#         'log_scale': False,
#         'cuts': {'sc_nClusters': lambda x: x > 0}  # Only superclusters with >0 clusters
#     },
#     {
#         'branches': ['sc_time', 'sc_energy'],
#         'name': 'sc_energy_vs_time',
#         'xlabel': 'SuperCluster Time [ns]',
#         'ylabel': 'SuperCluster Energy [MeV]',
#         'binning': [np.linspace(4, 27, 31), np.linspace(0, 1, 21)],
#         'log_scale': False,
#         'cuts': {'sc_nClusters': lambda x: x > 0}  # Only superclusters with >0 clusters
#     },

#     {
#         'branches': ['sc_nClusters', 'sc_energy'],
#         'name': 'sc_energy_vs_nClusters',
#         'xlabel': 'Number of Clusters',
#         'ylabel': 'SuperCluster Energy [MeV]',
#         'binning': [np.arange(-0.5, 3.5, 1), np.linspace(0, 1, 21)],
#         'log_scale': True,
#         'cuts': {'sc_nClusters': lambda x: x > 0}  # Only superclusters with >0 clusters
#     }


#     {
#         'branches': ['sc_nClusters', 'sd_primary_et'],
#         'name': 'sc_primaryPt_vs_nClusters',
#         'xlabel': 'Number of Clusters',
#         'ylabel': 'Primary particle pT [GeV]',
#         'binning': [np.arange(-0.5, 3.5, 1), np.linspace(0, 11, 21)],
#         'cuts': {'sc_nClusters': lambda x: x > 0}  # Only superclusters with >0 clusters
#     },

#     {
#         'branches': ['sc_nClusters', 'sd_primary_energy'],
#         'name': 'sc_primaryEnergy_vs_nClusters',
#         'xlabel': 'Number of Clusters',
#         'ylabel': 'Primary particle energy [GeV]',
#         'binning': [np.arange(-0.5, 3.5, 1), np.linspace(0, 40, 21)],
#         'cuts': {'sc_nClusters': lambda x: x > 0}  # Only superclusters with >0 clusters
#     },

#     {
#         'branches': ['sc_energy', 'sd_primary_energy'],
#         'name': 'sc_energy_vs_nClusters',
#         'xlabel': 'SC energy [MeV]',
#         'ylabel': 'Primary particle energy [GeV]',
#         'binning': [np.linspace(0, 5, 50), np.linspace(0, 40, 21)],
#         'cuts': {'sc_nClusters': lambda x: x > 0}  # Only superclusters with >0 clusters
#     },

# ]
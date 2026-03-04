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
    # NOTE: will throw StopIteration error if name not found
    next(filter(lambda x: x.get("name", None) == name, PLOT_DEFINITIONS_HISTORY))['binning'] = binning

PLOT_DEFINITIONS = PLOT_DEFINITIONS_HISTORY

# Example overlay with different cuts for the overlapped plot. Here we overlay
# the event-level 'simmc_time' but apply a tighter cut (nClusters > 2) to the
# overlaid dataset using an inline definition (dict) so we don't need a separate
# named PLOT_DEFINITION.
PLOT_DEFINITIONS.append({
    'branches': ['simmc_time_perCluster'],
    'name': 'simmc_time_perCluster_with_cut_variants',
    'xlabel': 'Cluster time [ns]',
    'ylabel': 'Clusters / MergedClusters',
    'binning': np.linspace(0, 30, 51),
    'type_breakdown': True,
    'cuts': {'simmc_nClusters': lambda x: x > 0},
    'overlap': [
        # named reference (same as before)
        'simmc_time',
        # inline dict: same branch but different cut
        {
            'branches': ['simmc_time'],
            'name': 'simmc_time_nClusters_gt2',
            'binning': np.linspace(0, 30, 51),
            'cuts': {'simmc_nClusters': lambda x: x > 2}
        }
    ],
})

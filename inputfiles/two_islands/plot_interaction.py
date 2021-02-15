"""
    PLOTS: - 'Relaxed magnetization angle' AS FUNCTION OF 'Time'
    Variable 'fancy': if True, angles are plotted in a continuous manner past [-180°, 180°].
                      if False, angles are constrained to [-180°, 180°].
    Variable 'groupBy': can be used to group by occurences of one table column (doesn't work entirely).
"""
import math
import matplotlib
import os
import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
from matplotlib import cm

font = {'size':16}
matplotlib.rc('font', **font)


def read_mumax3_table(filename):
    """Puts the mumax3 output table in a pandas dataframe"""
    
    table = pd.read_csv(filename, sep='\t')
    table.columns = ' '.join(table.columns).split()[1::2]
    
    return table

if __name__ == "__main__":
    # inFileName = 'two_islands_interaction.out/tableInt_a0Pi,0Pi_d128_r0.66,0.66_cell4nm.txt'
    # inFileName = 'two_islands_interaction.out/tableInt_a0Pi,0Pi_d128_r0.49,0.49_cell4nm.txt'
    inFileName = 'two_islands_interaction.out/tableInt_a0Pi,0Pi_d128_r0.81,0.81_cell4nm.txt'
    # inFileName = 'two_islands_interaction_noCustomField.out/table.txt'
    outDir = 'Figures/EnergyLandscape'
    if not os.path.exists(outDir):
        os.makedirs(outDir)
    outFileName = os.path.join(outDir, os.path.splitext(os.path.basename(inFileName).split('table')[-1])[0]) + '.pdf'

    table = read_mumax3_table(inFileName)

    fig = plt.figure(figsize=(7.0, 5.0))
    legend = []
    USE_ELECTRONVOLT = True
    USE_ABSOLUTE_VALUE = False
    WRAP_EDGES = True
    INTERPOLATE_PIXELS = True
    Energy = []
    if WRAP_EDGES:
        lim_1 = [0, 360]
        lim_2 = lim_1
    else:
        lim_1 = np.array([min(table["island1_angle"]), max(table["island1_angle"])])*180/np.pi
        lim_2 = np.array([min(table["island2_angle"]), max(table["island2_angle"])])*180/np.pi

    for subset in table.groupby("island1_angle", sort=True):
        island1_angle = subset[0]
        subtable = subset[1]
        if "E_Zeeman" in subtable.columns:
            diff = subtable["E_total"] - subtable["E_Zeeman"]
        else:
            diff = subtable["E_total"] - subtable["E_custom"]
        Energy.append((diff/1.602e-19 if USE_ELECTRONVOLT else diff))
    Energy = Energy - np.min(Energy)
    if WRAP_EDGES:
        Energy = np.insert(Energy, Energy.shape[0], values=Energy[0,:], axis=0)
        Energy = np.insert(Energy, Energy.shape[1], values=Energy[:,0], axis=1)
        

    ax = fig.add_subplot(111)
    interpolation = 'gaussian' if INTERPOLATE_PIXELS else 'nearest'
    im = ax.imshow(Energy[::-1], extent=[lim_1[0], lim_1[1], lim_2[0], lim_2[1]], interpolation=interpolation, cmap=cm.get_cmap('inferno'))
    # ax.set_aspect('auto')
    plt.xlabel("Island 2 magnetization angle [°]")
    plt.ylabel("Island 1 magnetization angle [°]")
    cbar = fig.colorbar(im)
    cbar.set_label('Energy [%s]' % ('eV' if USE_ELECTRONVOLT else 'J'), rotation=270, labelpad=25)

    angle_interval = 90
    plt.xticks(np.arange(lim_2[0], lim_2[1]+1e-3, angle_interval))
    plt.yticks(np.arange(lim_1[0], lim_1[1]+1e-3, angle_interval))

    plt.gcf().subplots_adjust(bottom=0.15, left=0.15, right=0.95)
    plt.gcf().tight_layout()
    plt.savefig(outFileName)

    plt.show()
    
"""
    PLOTS: - 'Height of energy landscape' AS FUNCTION OF 'Distance between center of islands'
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
    inFileName = 'two_islands_energyBarrier_distance.out/tableBarrierDistance_100-400_s100_r0.49_4nm.txt'
    # inFileName = 'two_islands_energyBarrier_distance.out/tableBarrierDistance_100-300_s100_r0.81_4nm.txt'
    outDir = 'Figures/EnergyBarrier'
    if not os.path.exists(outDir):
        os.makedirs(outDir)
    outFileName = os.path.join(outDir, os.path.splitext(os.path.basename(inFileName).split('table')[-1])[0]) + '.pdf'

    table = read_mumax3_table(inFileName)

    fig = plt.figure(figsize=(8.0, 5.0))
    legend = []
    USE_ELECTRONVOLT = True
    USE_ABSOLUTE_VALUE = False
    PLOT_AS_SCATTER = True
    GROUP_BY = "Cell_size"
    for subset in table.groupby(GROUP_BY, sort=False): # First loop can be ignored, perhaps can be used later
        size = subset[0]
        subtable = subset[1]
        E_barrier = []
        distances = []
        for subsubset in subtable.groupby("islands_distance"):
            distance = subsubset[0]
            subsubtable = subsubset[1]
            if "E_Zeeman" in subtable.columns:
                E_demag = subsubtable["E_total"] - subsubtable["E_Zeeman"]
            else:
                E_demag = subsubtable["E_total"] - subsubtable["E_custom"]
            if USE_ABSOLUTE_VALUE:
                diff = max(E_demag) - min(E_demag)
            elif len(E_demag) == 2:
                diff = E_demag.iloc[1] - E_demag.iloc[0] # 90° minus 0° value
            else:
                continue
            E_barrier.append(diff/1.602e-19 if USE_ELECTRONVOLT else diff)
            distances.append(distance)
        distances = np.array(distances)
        if GROUP_BY == "Size":
            legend.append('%s nm' % size)
        elif GROUP_BY == "Cell_size":
            legend.append('%s nm' % (size*1e9))
        if PLOT_AS_SCATTER:
            plt.scatter(distances*1e9, E_barrier)
        else:
            plt.plot(distances*1e9, E_barrier)

    plt.xlabel("Islands center distance [nm]")
    plt.ylabel("Energy landscape height [eV]")
    if len(legend) > 1:
        plt.legend(legend)
    plt.gcf().tight_layout()
    # plt.savefig(outFileName)

    plt.show()
    
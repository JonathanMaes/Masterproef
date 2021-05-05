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


def plot(inFileName, show=True, save=False, outName=''):
    table = read_mumax3_table(inFileName)

    USE_ELECTRONVOLT = True
    USE_ABSOLUTE_VALUE = False
    GROUP_BY = "size"
    highest_energy = 0
    fig = plt.figure(figsize=(8.0, 5.0))
    ax = fig.add_subplot(111)
    for subset in table.groupby(GROUP_BY, sort=False): # First loop can be ignored, perhaps can be used later
        size = subset[0]
        subtable = subset[1]
        E_barrier = []
        distancesNormalized = []
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
            distancesNormalized.append(distance/subsubtable["size"].iloc[0])
        distancesNormalized = np.array(distancesNormalized)
        highest_energy = max(max(E_barrier), highest_energy)

        if GROUP_BY == "size":
            label = '$L$ = %s nm' % int(size*1e9)
        elif GROUP_BY == "Cell_size":
            label = '%s nm' % (size*1e9)
            
        ax.plot(distancesNormalized, E_barrier, 'o', label=label)

        #### Fit to 1/dist³
        lastValue = E_barrier[-1]
        lastDistance = distancesNormalized[-1]
        fittedEnergies = [lastValue*(dist/lastDistance)**-3 for dist in distancesNormalized]
        ax.plot(distancesNormalized, fittedEnergies, '--', color='black', linewidth=1, label=None)

    plt.xlabel(r"Normalized distance $d/L$")
    plt.ylabel("Energy landscape height [eV]")
    plt.ylim([highest_energy*(-0.01), highest_energy*1.1])
    plt.xlim([0.97, np.max(distancesNormalized)+0.03])
    plt.legend()
    
    if show:
        plt.gcf().tight_layout()
        if save:
            outDir = 'Figures/EnergyBarrier'
            if not os.path.exists(outDir):
                os.makedirs(outDir)
            if outName:
                outFileName = os.path.join(outDir, outName) + '.pdf'
            else:
                outFileName = os.path.join(outDir, os.path.splitext(os.path.basename(inFileName).split('table')[-1])[0]) + '.pdf'
            plt.savefig(outFileName)
        plt.show()


if __name__ == "__main__":
    plot('two_islands_energyBarrier_distance.out/tableBarrierDistance_s50,100_r0.49_dist1-4L_2,4nm.txt', save=True, outName='dist1-4L_r0.49_s100&50_4nm&2nm')

    # plot('two_islands_energyBarrier_distance.out/tableBarrierDistance_s100_r0.81_dist100-300_4nm.txt')
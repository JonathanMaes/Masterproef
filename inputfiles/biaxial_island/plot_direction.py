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

font = {'size':16}
matplotlib.rc('font', **font)


def read_mumax3_table(filename):
    """Puts the mumax3 output table in a pandas dataframe"""
    
    table = pd.read_csv(filename, sep='\t')
    table.columns = ' '.join(table.columns).split()[1::2]
    
    return table

if __name__ == "__main__":
    # inFileName = 'biaxial_island_switching_plus.out/table_65x100_300K_alpha0.1_1µs.txt'
    # inFileName = 'biaxial_island_switching_plus.out/table_65x100_300K_alpha0.01_1µs.txt'
    inFileName = 'biaxial_island_switching_plus.out/table_100x100_300K_alpha0.01_1µs.txt'
    outDir = 'Figures/Switching'
    if not os.path.exists(outDir):
        os.makedirs(outDir)
    outFileName = os.path.join(outDir, os.path.splitext(os.path.basename(inFileName).split('table_')[-1])[0]) + '.pdf'

    table = read_mumax3_table(inFileName)

    fancy = True # If fancy: follow angles past [-180°, 180°]
    groupBy = "" # If groupBy is something: multiple plots grouped by that property
    if groupBy:
        subsets = [subset[1] for subset in table.groupby(groupBy)]
        legend = ['%s: %s' % (groupBy, subset[0]) for subset in table.groupby(groupBy)]
    else:
        subsets = [table]

    fig = plt.figure(figsize=(8.0, 5.0))
    for subtable in subsets:
        angles = np.arctan2(subtable["my"], subtable["mx"])*180/math.pi
        if not fancy:
            plt.plot(subtable["t"]*1e9, angles)
        else:
            previousAngles = np.array(angles)[:-1]
            nextAngles = np.array(angles)[1:]
            offsets = (np.abs(previousAngles - nextAngles) > 180)*((previousAngles > nextAngles)*2 - 1)*360
            offset = np.cumsum(offsets)
            fancyAngles = np.append([angles[0]], angles[1:] + offset)
            plt.plot(subtable["t"]*1e9, fancyAngles)

    if groupBy:
        plt.legend(legend)

    plt.xlabel(r'$t$ [ns]')
    plt.ylabel(r'angle [°]')
    plt.gcf().tight_layout()
    plt.savefig(outFileName)

    plt.show()
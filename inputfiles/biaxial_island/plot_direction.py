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
    # inFileName = 'biaxial_island_switching_plus.out/table_65x100_300K_alpha0.1_1µs_4nm.txt'
    inFileName = 'biaxial_island_switching_plus.out/table_65x100_300K_alpha0.01_1µs_4nm.txt'
    # inFileName = 'biaxial_island_switching_plus.out/table_65x100_350K_alpha0.01_1µs_4nm.txt'
    # inFileName = 'biaxial_island_switching_plus.out/table_65x100_273K_alpha0.01_1µs_4nm.txt'
    # inFileName = 'biaxial_island_switching_plus.out/table_100x100_350K_alpha0.01_1µs_4nm.txt'
    # inFileName = 'biaxial_island_switching_plus.out/table_49x100_300K_alpha0.01_100ns_2nm.txt'
    # inFileName = 'biaxial_island_switching_plus.out/table_65x100_350K_alpha0.01_0.5µs_2nm.txt'
    # inFileName = 'biaxial_island_switching_plus.out/table_65x100_350K_alpha0.01_1µs_3.125nm.txt'
    # inFileName = 'biaxial_island_switching_extfield.out/table_50x45_ext0.00015_ext3Pi8_1.25µs_3.125nm.txt'
    # inFileName = 'biaxial_island_switching_extfield.out/table_50x45_ext0_10ns_3.125nm.txt'
    # inFileName = 'biaxial_island_switching_extfield.out/table.txt'
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
    angleRange = [0, 0]
    timeRange = [0, 0]
    for subtable in subsets:
        angles = np.arctan2(subtable["my"], subtable["mx"])*180/math.pi
        timeRange[0] = min(timeRange[0], np.min(subtable["t"])*1e9)
        timeRange[1] = max(timeRange[1], np.max(subtable["t"])*1e9)
        if not fancy:
            plt.plot(subtable["t"]*1e9, angles)
        else:
            previousAngles = np.array(angles)[:-1]
            nextAngles = np.array(angles)[1:]
            offsets = (np.abs(previousAngles - nextAngles) > 180)*((previousAngles > nextAngles)*2 - 1)*360
            offset = np.cumsum(offsets)
            fancyAngles = np.append([angles[0]], angles[1:] + offset)
            angleRange[0] = min(angleRange[0], np.min(fancyAngles))
            angleRange[1] = max(angleRange[1], np.max(fancyAngles))
            plt.plot(subtable["t"]*1e9, fancyAngles)

    if groupBy:
        plt.legend(legend)

    plt.grid(axis='y', color='grey', linestyle=':', linewidth=1)
    plt.xlabel(r'$t$ [ns]')
    plt.ylabel(r'angle [°]')
    plt.xlim(timeRange[0], timeRange[1])
    angle_interval = 90
    plt.yticks(np.arange(angleRange[0]//angle_interval*angle_interval, angleRange[1]//angle_interval*angle_interval+angle_interval+1, angle_interval))
    plt.gcf().tight_layout()
    plt.savefig(outFileName)

    plt.show()
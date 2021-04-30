"""
    PLOTS: - 'Energy barrier' AS FUNCTION OF 'Roundness' FOR DIFFERENT 'Ellipse size'
"""
import matplotlib
import os
import numpy as np
import matplotlib.pyplot as plt
import pandas as pd

font = {'size':16}
matplotlib.rc('font', **font)
matplotlib.rcParams['axes.prop_cycle'] = matplotlib.cycler(color=["b", "y", "r", "lightgreen"]) 


def read_mumax3_table(filename):
    """Puts the mumax3 output table in a pandas dataframe"""
    
    table = pd.read_csv(filename, sep='\t')
    table.columns = ' '.join(table.columns).split()[1::2]
    
    return table

if __name__ == "__main__":
    inFileName = 'biaxial_island_shape.out/tablePlus_100_0.1-1_aPi4_B0.001.txt' #! set GROUP_BY to "Cell_size" for this one
    inFileName = 'biaxial_island_shape.out/tablePlus_100_0.1-1.txt' #! set GROUP_BY to "Cell_size" for this one
    # inFileName = 'biaxial_island_shape.out/tablePlus_32-128_0.1-1_aPi4_B0.001_cell1nm.txt' #! set GROUP_BY to "Size" for this one
    # inFileName = 'biaxial_island_shape.out/tablePlus_100_0.1-1_aPi4_B0.001_cell3.125nm.txt'
    inFileName = 'biaxial_island_shape.out/tablePlus_100_0.1-1_aPi128_B0.01_cell4nm.txt' # Those with aPi128 include the second-order anisotropy, at the cost of using a higher field.
    inFileName = 'biaxial_island_shape.out/table.txt'
    outDir = 'Figures/Barrier'
    if not os.path.exists(outDir):
        os.makedirs(outDir)
    outFileName = os.path.join(outDir, os.path.splitext(os.path.basename(inFileName).split('table')[-1])[0]) + '.pdf'

    shape = read_mumax3_table(inFileName)

    fig = plt.figure(figsize=(8.0, 5.0))
    legend = []
    USE_ELECTRONVOLT = True
    GROUP_BY = "Size"
    # GROUP_BY = "Cell_size"
    USE_ABSOLUTE_VALUE = False
    PLOT_AS_SCATTER = False
    for subset in shape.groupby(GROUP_BY, sort=False):
        size = subset[0]
        subtable = subset[1]
        E_barrier = []
        roundnesses = []
        for subsubset in subtable.groupby("Roundness"):
            roundness = subsubset[0]
            subsubtable = subsubset[1]
            E_demag = subsubtable["E_total"]-subsubtable["E_Zeeman"]
            if USE_ABSOLUTE_VALUE:
                diff = max(E_demag) - min(E_demag)
            elif len(E_demag) == 2:
                diff = E_demag.iloc[1] - E_demag.iloc[0]
            else:
                allowed_energies = E_demag
                allowed_energies = [w for i, w in enumerate(E_demag) if subsubtable["my"].iloc[i] >= -1e-6 and subsubtable["mx"].iloc[i] >= -1e-6]
                # Find index closest to 0°
                index_0 = (np.abs(subsubtable["my"])).idxmin()
                # Find index closest to 45°
                index_45 = (np.abs(np.abs(subsubtable["my"]) - np.sqrt(2)/2)).idxmin()
                # print(E_demag, index_0, index_45)
                # Determine if energy barrier is positive (easy axis 0°) or negative (easy axis 45°)
                sign = 1 if E_demag[index_0] < E_demag[index_45] else -1
                diff = (max(allowed_energies) - min(allowed_energies))*sign
                print(diff)
            # print("(%.2f x %.2f nm) Delta E: %.2e J = %.3f eV" % (size, size*roundness, diff, diff/1.602e-19))
            E_barrier.append(diff/1.602e-19 if USE_ELECTRONVOLT else diff)
            roundnesses.append(roundness)
        if GROUP_BY == "Size":
            legend.append('%s nm' % size)
        elif GROUP_BY == "Cell_size":
            legend.append('%s nm' % (size*1e9))
        if PLOT_AS_SCATTER:
            plt.scatter(roundnesses, E_barrier)
        else:
            plt.plot(roundnesses, E_barrier)

    plt.grid(color='grey', linestyle=':', linewidth=1)
    plt.xlabel(r'Roundness $\rho$')
    plt.ylabel(r'Energy barrier [%s]' % ('eV' if USE_ELECTRONVOLT else 'J'))
    plt.legend(legend)
    plt.gcf().tight_layout()
    # plt.savefig(outFileName)
    plt.show()
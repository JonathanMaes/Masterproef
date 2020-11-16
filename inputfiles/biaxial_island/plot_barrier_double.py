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


def read_mumax3_table(filename):
    """Puts the mumax3 output table in a pandas dataframe"""
    
    table = pd.read_csv(filename, sep='\t')
    table.columns = ' '.join(table.columns).split()[1::2]
    
    return table

if __name__ == "__main__":
    inFileName = 'biaxial_island_shape.out/tablePlus_100_0.1-1_aPi4_B0.001_cell2nm.txt'
    # inFileName = 'biaxial_island_shape.out/tablePlus_16-128_0.1-1_aPi4_B0.001_cell2nm.txt'
    outDir = 'Figures/Barrier'
    if not os.path.exists(outDir):
        os.makedirs(outDir)
    outFileName = os.path.join(outDir, os.path.splitext(os.path.basename(inFileName).split('table')[-1])[0]) + '.pdf'

    shape = read_mumax3_table(inFileName)

    fig = plt.figure(figsize=(8.0, 5.0))
    legend = []
    USE_ELECTRONVOLT = True
    USE_ABSOLUTE_VALUE = False
    for subset in shape.groupby("Size", sort=False):
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
                continue
            # print("(%.2f x %.2f nm) Delta E: %.2e J = %.3f eV" % (size, size*roundness, diff, diff/1.602e-19))
            E_barrier.append(diff/1.602e-19 if USE_ELECTRONVOLT else diff)
            roundnesses.append(roundness)
        legend.append('%d nm' % size)
        plt.plot(roundnesses, E_barrier)

    plt.grid(color='grey', linestyle=':', linewidth=1)
    plt.xlabel(r'Roundness')
    plt.ylabel(r'Energy barrier [%s]' % ('eV' if USE_ELECTRONVOLT else 'J'))
    plt.legend(legend)
    plt.gcf().tight_layout()
    plt.savefig(outFileName)
    plt.show()
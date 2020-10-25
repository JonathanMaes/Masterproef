"""
    PLOTS: - 'Energy barrier' AS FUNCTION OF 'Roundness' FOR DIFFERENT 'Ellipse size'
"""
import numpy as np
import matplotlib.pyplot as plt
import pandas as pd


def read_mumax3_table(filename):
    """Puts the mumax3 output table in a pandas dataframe"""
    
    table = pd.read_csv(filename, sep='\t')
    table.columns = ' '.join(table.columns).split()[1::2]
    
    return table

if __name__ == "__main__":
    # Shape anisotropy: Plus
    shape = read_mumax3_table('biaxial_island_shape.out/tablePlus_10-100_0.1-1_aPi4_B0.001_l5_0.05.txt')
    # shape = read_mumax3_table('biaxial_island_shape.out/tablePlus_128-16_0.1-1_aPi4_B0.001.txt')
    legend = []
    for subset in shape.groupby("Size"):
        size = subset[0]
        subtable = subset[1]
        E_barrier = []
        roundnesses = []
        for subsubset in subtable.groupby("Roundness"):
            roundness = subsubset[0]
            subsubtable = subsubset[1]
            E_demag = subsubtable["E_total"]-subsubtable["E_Zeeman"]
            diff = max(E_demag) - min(E_demag)
            print("(%.2f x %.2f nm) Delta E: %.2e J = %.3f eV" % (size, size*roundness, diff, diff/1.602e-19))
            E_barrier.append(diff)
            roundnesses.append(roundness)
        legend.append('%d nm' % size)
        plt.plot(roundnesses, E_barrier)
    plt.title("Double ellipse geometry")
    plt.legend(legend)
    plt.xlabel("Roundness [nm]")
    plt.ylabel("Energy barrier [J]")
    plt.show()
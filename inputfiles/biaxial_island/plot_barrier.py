"""
    PLOTS: - 'Energy' AS FUNCTION OF 'External magnetic field angle' FOR DIFFERENT 'Ellipse width'
           - 'Energy barrier' AS FUNCTION OF 'Ellipse width'
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
    # Cubic anisotropy
    cubic = read_mumax3_table('biaxial_island_cubic.out/table.txt')
    legend = ['Cubic']
    plt.plot(cubic["Angle"], cubic["E_total"])

    # Shape anisotropy: Plus
    # shape = read_mumax3_table('biaxial_island_shape.out/tablePlus_10-100_aPi4_B0.001_l1.txt')
    shape = read_mumax3_table('biaxial_island_shape.out/tablePlus_10-100_a0-2Pi_Pi4_B0.001_l5.txt')
    E_barrier = []
    sizes = []
    for subset in shape.groupby("Length"):
        size = subset[0]
        subtable = subset[1]
        E_demag = subtable["E_total"]-subtable["E_Zeeman"]
        diff = max(E_demag) - min(E_demag)
        print("(%.2f nm) Delta E: %.2e J = %.3f eV" % (size, diff, diff/1.602e-19))
        E_barrier.append(diff)
        sizes.append(size)
        legend.append('Shape%d' % size)
        plt.plot(subtable["Angle"], E_demag)
    plt.legend(legend)
    plt.title("Double ellipse geometry")
    plt.xlabel("External magnetic field angle [°]")
    plt.ylabel("Energy [J]")
    plt.show()

    # Barrier as function of size
    plt.plot(sizes, E_barrier)
    plt.title("Double ellipse geometry")
    plt.xlabel("Ellipse width [nm]")
    plt.ylabel("Energy barrier [J]")
    plt.show()
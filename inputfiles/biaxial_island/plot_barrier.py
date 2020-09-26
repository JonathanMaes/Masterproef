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
    shape = read_mumax3_table('biaxial_island_shape.out/tablePlus_70-100_0.03step.txt')
    for subset in shape.groupby("Length"):
        size = subset[0]
        subtable = subset[1]
        E_demag = subtable["E_total"]-subtable["E_Zeeman"]
        diff = max(E_demag) - min(E_demag)
        print("(%.2f nm) Delta E: %.2e J = %.3f eV" % (size, diff, diff/1.602e-19))
        legend.append('Shape%d' % size)
        plt.plot(subtable["Angle"], E_demag)
    plt.legend(legend)

    # Show plot
    plt.show()
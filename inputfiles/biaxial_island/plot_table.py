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

    # Shape anisotropy
    shape = read_mumax3_table('biaxial_island_shape.out/table.txt')
    for subset in shape.groupby("Length"):
        size = subset[0]
        subtable = subset[1]
        legend.append('Shape%d' % size)
        plt.plot(subtable["Angle"], subtable["E_total"])
    plt.legend(legend)

    # Show plot
    plt.show()
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
    shape = read_mumax3_table('biaxial_island_shape_field.out/tablePlus_65_T100-1.75_a0.1.txt')
    legend = []
    for subset in shape.groupby("Field"):
        field = subset[0]
        subtable = subset[1]
        E_demag = subtable["E_total"]-subtable["E_Zeeman"]
        diff = max(E_demag) - min(E_demag)
        print("(%.2f nm) Delta E: %.2e J = %.3f eV" % (field, diff, diff/1.602e-19))
        legend.append('%.2f T' % field)
        plt.plot(subtable["Angle"], E_demag)
    plt.legend(legend)

    # Show plot
    plt.show()
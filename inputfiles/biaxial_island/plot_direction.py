import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
import math


def read_mumax3_table(filename):
    """Puts the mumax3 output table in a pandas dataframe"""
    
    table = pd.read_csv(filename, sep='\t')
    table.columns = ' '.join(table.columns).split()[1::2]
    
    return table

if __name__ == "__main__":
    table = read_mumax3_table('biaxial_island_switching_plus_remote.out/table_l65_1µs.txt')
    angles = np.arctan2(table["my"], table["mx"])*180/math.pi

    fancy = True # If fancy: follow angles past [-180°, 180°]
    if not fancy:
        plt.plot(table["t"], angles)
    else:
        previousAngles = np.array(angles)[:-1]
        nextAngles = np.array(angles)[1:]
        offsets = (np.abs(previousAngles - nextAngles) > 180)*((previousAngles > nextAngles)*2 - 1)*360
        offset = np.cumsum(offsets)
        fancyAngles = np.append([angles[0]], angles[1:] + offset)
        plt.plot(table["t"], fancyAngles)

    plt.xlabel("t [s]")
    plt.ylabel("angle [°]")

    plt.show()
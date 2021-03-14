"""
    PRINTS
    Values of the energy minima for each of four possible input directions,
    and the corresponding relaxed magnetization angles of the other islands, in order of their region numbers.
"""


import math
import matplotlib
import os
import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
from matplotlib import cm

font = {'size':16}
matplotlib.rc('font', **font)


def read_mumax3_table(filename):
    """Puts the mumax3 output table in a pandas dataframe"""
    
    table = pd.read_csv(filename, sep='\t')
    table.columns = ' '.join(table.columns).split()[1::2]
    
    return table

def convert_ovf(folder):
    os.system('mumax3-convert -png -arrows 8 %s/*.ovf' % folder)

def anglediff(a1, a2):
    return (a1 - a2 + 180 + 360) % 360 - 180

def process_data(filename, input_island=1):
    table = read_mumax3_table(filename)
    num_islands = 0
    for name in table.columns:
        if name.startswith('m.region') and name.endswith('x'):
            num_islands += 1

    island_magAngles = []
    for i in range(num_islands):
        island_magAngles.append(np.arctan2(table['m.region%dy' % (i+1)], table['m.region%dx' % (i+1)])*180/np.pi)
    island_magAngles = np.array(island_magAngles).transpose()
    energies = np.array(table["E_total"] - table["E_Zeeman"])

    # sorted_island_magAngles = island_magAngles[energies.argsort(),:]
    # sorted_energies = energies[energies.argsort()]

    a1_minEnergies = []
    a1_minEnergyMagAngles = []
    for i, subset in enumerate(table.groupby("a%d" % input_island, sort=False)):
        a1 = subset[0]
        # subtable = subset[1]

        input_island_angles = island_magAngles[:,input_island-1].transpose()

        indices_close_to_a1 = np.where(np.abs(anglediff(input_island_angles, a1*180/np.pi)) < 22.5)[0]
        if len(indices_close_to_a1) == 0:
            print('No stable configurations for input at %.2f rad' % a1)

        a1_energies = energies[indices_close_to_a1]
        a1_magAngles = island_magAngles[indices_close_to_a1]
        
        sorted_a1_magAngles = a1_magAngles[a1_energies.argsort(),:]
        sorted_a1_energies = a1_energies[a1_energies.argsort()]

        degeneracy = np.count_nonzero(sorted_a1_energies == sorted_a1_energies[0])

        a1_minEnergies.append(sorted_a1_energies[0])
        a1_minEnergyMagAngles.append(sorted_a1_magAngles[:degeneracy])

    print('Lowest energies: \n', a1_minEnergies)
    print('\nStable angles: \n', np.array(a1_minEnergyMagAngles))






if __name__ == "__main__":
    convert_ovf('many_islands_interaction.out')
    print('#'*80)

    process_data('many_islands_interaction.out/table.txt', input_island=1)
    # process_data('attempts/table000001.txt', input_island=1)

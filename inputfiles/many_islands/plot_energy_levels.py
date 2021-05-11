"""
    This Python program uses the ++ half adder, and modifies the parameters of the fixed island.
    Since this requires determining if the configuration is a half adder, for many different simulations,
    the program will call check_if_halfadder on subsets of a huge table.txt.
"""

import matplotlib
import matplotlib.pyplot as plt
import numpy as np
import os
import pandas as pd
import shutil
from math import pi
from matplotlib import cm
from matplotlib.colors import ListedColormap, LinearSegmentedColormap
from matplotlib.patches import Rectangle
font = {'size':16}
matplotlib.rc('font', **font)

from read_output_manyislands import check_if_halfadder, read_mumax3_table, process_data, collect_orientation, round_90
from generate_mumax3_inputfile import generate_mumax3_inputfile, Island

# Define some filenames
# attemptcode = 'test'
attemptcode = '000006'
outDir = 'Results/Halfadder%s' % attemptcode
if not os.path.exists(outDir):
    os.makedirs(outDir)


def replace_with_dict(ar, dic):
    ''' Replaces the values of the numpy array <ar> according to the rules in the dictionary <dic>. '''
    k = np.array(list(dic.keys()))
    v = np.array(list(dic.values()))
    sidx = k.argsort()
    return v[sidx[np.searchsorted(k,ar,sorter=sidx)]]

def plot_energy_levels(halfadderfile):
    """
        @param halfadderfile [str]: The mumax table.txt containing all possible states of the half adder.
    """
    
    #### Read table and determine all relevant half adder quantities
    table = read_mumax3_table(halfadderfile)

    halfadder = check_if_halfadder(table, verbose=False)[0]
    mag_angles, energies, geom_angles = process_data(table)

    input_island = halfadder[0]
    output_island = 2 if input_island == 1 else 1
    halfadder_comb = halfadder[3]

    collected_energies = []
    collected_magangles = []
    largest_collection = 0 # This is the largest number of energy levels
    for i in range(4):
        new_mag_angles, new_energies = collect_orientation(input_island, geom_angles[input_island-1]+i*90, mag_angles, energies, geom_angles)
        collected_energies.append(new_energies)
        collected_magangles.append(new_mag_angles[:,output_island-1])
        largest_collection = max(largest_collection, len(collected_energies))
    collected_energies_matrix = np.array([list(l) + [np.infty]*(largest_collection - len(l)) for l in collected_energies])

    highest_groundstate = np.max(collected_energies_matrix[:,0])
    lowest_excitedstate = np.min(collected_energies_matrix[:,1])
    print(highest_groundstate, lowest_excitedstate)

    print(collected_energies, collected_magangles)
    energyheight = np.max([np.max(collected_energies[i]) for i in range(4)])
    print(energyheight)
    print(np.max(collected_energies_matrix))

    #### Plots    
    fig = plt.figure(figsize=(10.0, 5.0))
    ax = fig.add_subplot(111)

    PLOTVAR_linewidth = 0.8
    PLOTVAR_bboxprops = dict(boxstyle='round', facecolor='wheat', alpha=1)
    PLOTVAR_minallowedenergydiff_fraction = 0.05

    ax.axhline(highest_groundstate, color='blue', linestyle='--', linewidth=1)
    ax.axhline(lowest_excitedstate, color='black', linestyle='--', linewidth=1)
    
    resulting_comb = []
    for i in range(4):
        current_energies = collected_energies[i]
        current_magangles = collected_magangles[i]
        current_magangles_clean = round_90(current_magangles, delta=geom_angles[output_island-1]) % 360
        print(current_magangles)
        print(current_magangles_clean)
        if len(current_energies) > 1:
            minenergydiff = min(current_energies[1:]-current_energies[:-1])
        else:
            minenergydiff = np.inf
        print(minenergydiff)
        for j, energylevel in enumerate(current_energies):
            ax.plot([i+(1-PLOTVAR_linewidth)/2, i+1-(1-PLOTVAR_linewidth)/2], [energylevel]*2, color='black' if j != 0 else 'blue', linewidth=3)
            text = '%d°'% current_magangles_clean[j]
            if j == 0:
                lowestenergy_comb = halfadder_comb[int(current_magangles_clean[j]/90)]
                # text += '=%d' % lowestenergy_comb
                resulting_comb.append(lowestenergy_comb)
            textpos_x = i+0.5
            if minenergydiff < energyheight*PLOTVAR_minallowedenergydiff_fraction: # If energy levels too close, then put the text off-center
                textpos_x = i+.2+0.2*j
            ax.text(textpos_x, energylevel, text, fontsize=10, transform=ax.transData, ha='center', va='center', bbox=PLOTVAR_bboxprops)
    # for line in contourlines:
    #     ax.plot(line[1], line[0], **kwargs)

    ax.tick_params(axis='both', which='major', length=8)
    plt.xticks([i+0.5 for i in range(4)], [('In %d°\n%d' + r' $\rightarrow$ %d') % (i*90, halfadder_comb[i], resulting_comb[i]) for i in range(4)])
    ax.set_ylabel('Energy [eV]')
    ax.spines['right'].set_color('none')
    ax.spines['top'].set_color('none')

    # plt.gcf().subplots_adjust(bottom=0.15, left=0.15, right=0.95)
    plt.gcf().tight_layout()
    plt.savefig(os.path.join(outDir, os.path.split(halfadderfile)[1].replace('.txt', '_energylevels.pdf')))
    plt.show()




if __name__ == "__main__":
    pass
    plot_energy_levels('Results/Halfadder000006/table(d170,s-60).txt')

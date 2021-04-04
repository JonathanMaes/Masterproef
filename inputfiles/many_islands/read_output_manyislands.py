"""
    PRINTS
    Values of the energy minima for each of four possible input directions,
    and the corresponding relaxed magnetization angles of the other islands, in order of their region numbers.
"""

import itertools
import math
import matplotlib
import matplotlib.pyplot as plt
import os
import re
import numpy as np
import pandas as pd
from matplotlib import cm

font = {'size':16}
matplotlib.rc('font', **font)
matplotlib.rcParams['axes.prop_cycle'] = matplotlib.cycler(color=["b", "y", "r", "lightgreen"]) 


def read_mumax3_table(filename):
    """Puts the mumax3 output table in a pandas dataframe"""
    
    table = pd.read_csv(filename, sep='\t')
    table.columns = ' '.join(table.columns).split()[1::2]
    
    return table

def convert_ovf(folder):
    os.system('mumax3-convert -png -arrows 8 %s/*.ovf' % folder)

def anglediff(a1, a2):
    return (a1 - a2 + 180 + 360) % 360 - 180

def round_90(a, delta=0):
    return np.round((a-delta)/90)*90+delta

def rad_to_deg(a):
    return a*180/math.pi
def deg_to_rad(a):
    return a*math.pi/180

def process_data(filename, sort_energies=True, normalize_energies=True):
    '''
        Reads mumax3 table at <filename> and gives for each attempt:
            - the magnetization angles of all islands
            - the corresponding energy
        and also gives the geometry angle of each island.

        @param filename: Relative path to mumax3 table.txt file.
        @param sort_energies (True): Whether the output energies and corresponding mag_angles should be sorted.
        @returns tuple(3):
            mag_angles [deg]: List of lists, with each sublist being one a{n} situation. Element {n} of a sublist is the relaxed m angle of island {n} for that situation.
            energies [eV]: List, with element n being the energy corresponding to mag_angles[n] relaxed configuration.
            geom_angles [deg]: The angle under which each geometry is rotated, in order.
            All these include fixed islands, and energies are not normalized to the minimum one.
    '''
    table = read_mumax3_table(filename)
    num_islands = 0
    for name in table.columns:
        # If column is of form "m.region{n}x"
        if re.compile(r"^m\.region\d+x$").match(name):
            num_islands += 1
    geom_angles = []
    for i in range(1,num_islands+1):
        if 'a%d' % i in table.columns:
            geom_angles.append(rad_to_deg(table['a%d'%i].iloc[0]) % 90)
        else: # Default for fixed islands because we dont care about these in the output anyway
            geom_angles.append(0)

    mag_angles = [] # List of lists, with each sublist being one a{n} situation. Element {n} of a sublist is the relaxed m angle of island {n} for that situation.
    for i in range(num_islands):
        mag_angles.append(rad_to_deg(np.arctan2(table['m.region%dy' % (i+1)], table['m.region%dx' % (i+1)])))
    mag_angles = np.array(mag_angles).transpose()
    energies = np.array(table["E_total"] - table["E_Zeeman"])

    if sort_energies:
        mag_angles = mag_angles[energies.argsort(),:]
        energies = energies[energies.argsort()]
    energies = energies/1.602e-19
    if normalize_energies:
        energies -= min(energies)

    return (mag_angles, energies, geom_angles)


def collect_orientation(island, angle, mag_angles, energies, geom_angles, margin_mag=45, margin_energy=0.01, sort_energies=True):
    """
        Given an island number <island>, and a desired direction <angle>, this function
        looks at only those entries which have island <island> at or near that angle.
        This can be used to find the minimal energy configuration if one were to use island
        number <island> as input bit and use as input the angle <angle>.

        @param island [int]: The island that is considered the 'input island'.
        @param angle [deg]: The angle under which the input island is forced.
        @params mag_angles, energies, geom_angles: Simply the output of process_data().
        @param margin_mag [deg] (22.5): Margin by which two magnetization angles are considered the same.
        @param margin_energy [eV] (0.01): Margin by which two energies are considered the same.
        @param sort_energies [bool] (True): Whether the output energies and corresponding mag_angles should be sorted.
    """
    new_mag_angles = []
    new_energies = []

    # Only take those values for mag_angles[island-1] which are close to <angle>
    indices_close_to_a1 = np.where(np.abs(anglediff(mag_angles[:,island-1].transpose(), angle)) < margin_mag)[0]
    # if len(indices_close_to_a1) == 0:
    #     print('No stable configurations for input island %d at %.2f deg' % (island, angle))
    energies = energies[indices_close_to_a1]
    mag_angles = mag_angles[indices_close_to_a1]

    for i, entry in enumerate(mag_angles):
        # Search for a similar entry, if both mag_angles and energies are similar, then discard this one
        if np.any(np.all(np.abs(anglediff(entry, mag_angles[:i])) < margin_mag, axis=1)) and np.any(np.abs(energies[i]-energies[:i]) < margin_energy):
            continue # Similar entry has already been recorded, so dont record <i>
        else:
            new_mag_angles.append(entry)
            new_energies.append(energies[i])
    return (np.array(new_mag_angles), np.array(new_energies))

# TODO: Now write a function that does
# - Extra function to plot energies corresponding to a certain input for a certain island
# - Perhaps with nice GUI but that might be a lot of work

def get_lowest_energies(filename, input_island, verbose=True):
    '''
        For the mumax3 table file at <filename>, this function lists the minimal energy configuration
        for each of 4 possible input bits if one considers island <input_island> to be the input island.
        Hence, the output is a list with 4 elements, each element corresponding to a different input bit for
        the input island. This element is then a list containing the relaxed magnetization angles of the
        lowest energy state for which island <input_island> is at an angle corresponding to that input bit.

        @param filename [str]: Relative path to mumax3 table.txt file.
        @param input_island [int]: The island that is considered the 'input island'.
        @param verbose [bool] (True): Whether to print the @returns content to the terminal or not.
        @returns: - None if no stable configurations exist for some angle of input_island, for example because
                    the island is fixed (then no stable configurations exist for any angle that is not the fixation direction)
                  - Otherwise, returns (minEnergyMagAngles, minEnergies), where
                    minEnergyMagAngles: A list with 4 elements, each corresponding to a different input bit for the input island.
                        The elements are themselves lists, containing the relaxed magnetization angles of the islands at the
                        minimal energy configuration. The corresponding energy is returned in <minEnergies>.
                    minEnergies: A list with 4 elements, corresponding to the energy of the respective configuration in <minEnergyMagAngles>.
    '''
    mag_angles, energies, geom_angles = process_data(filename)

    minEnergies = []
    minEnergyMagAngles = []
    for i in range(4):
        new_mag_angles, new_energies = collect_orientation(input_island, geom_angles[input_island-1]+i*90, mag_angles, energies, geom_angles)

        if len(new_energies) == 0: # No stable configurations for <input_island> at geom_angles[i]
            print('Island %d does not have a stable orientation around %d degrees.' % (input_island, geom_angles[input_island-1]+i*90))
            return (None, None)
        else:
            minEnergies.append(new_energies[0])
            minEnergyMagAngles.append(new_mag_angles[0])
            if abs(new_energies[0] - new_energies[1]) < 1e-4:
                print('Input island %d at %d deg is DEGENERATE!' % (input_island, geom_angles[input_island-1]+i*90))

    if not verbose:
        print('For island %d:' % input_island)
        print(' Lowest energies: \n %s' % minEnergies)
        print(' Stable angles: \n %s' % np.array(minEnergyMagAngles))
    return (np.array(minEnergyMagAngles), np.array(minEnergies))

def check_if_halfadder(filename, allow_rotation=False):
    """
        For a given mumax3 table file at <filename>, this function looks if there is a half adder
        possible with the relaxed magnetizations from the table file.

        @param filename [str]: Relative path to mumax3 table.txt file.
        @param allow_rotation [bool] (False): If true, the search does not care about the orientation of input
                                              and output with respect to each other, so this allows for 45°
                                              difference to be neglected. TODO!
        @returns: - False if no half adder was found.
                  - (i, j, angles_i, input_bits), where
                    i: Index of input island.
                    j: Index of output island.
                    angles_i: The four stable angles for the input island magnetization.
                    input_bits: The bit associated with each magnetization angle in angles_i.
    """
    _, _, geom_angles = process_data(filename)
    halfAdder = {0:0,1:1,2:1,3:2}
    found_halfadders = []
    for i in range(len(geom_angles)): # Input island <i>
        print('Checking if island %d can be an input...' % (i+1))
        minEnergyMagAngles, _ = get_lowest_energies(filename, i+1)
        if minEnergyMagAngles is None: # Then no stable configurations exist for some angle of input <i> (e.g. because island <i> is fixed)
            continue
        angles_i = anglediff(0, round_90(minEnergyMagAngles[:,i].T, delta=geom_angles[i]))
        if np.all(np.array(angles_i) == angles_i[0]): # Then island <i> is fixed
            continue
        for j in range(len(geom_angles)): # Output island <j>
            if i == j: # If same island, then skip
                continue
            if not allow_rotation:
                if geom_angles[i] != geom_angles[j]: # If not same orientation, then skip
                    continue
            angles_j = anglediff(0, round_90(minEnergyMagAngles[:,j].T, delta=geom_angles[j]))
            if np.all(np.array(angles_j) == angles_j[0]): # Then island <j> is fixed
                continue
            for perm in itertools.permutations([0,1,2,3]):
                d = {angles_i[n]:perm[n] for n in range(4)}
                input_bits = perm
                should_be_output_bits = list(map(halfAdder.get, input_bits))
                output_bits = list(map(d.get, angles_j))
                if should_be_output_bits == output_bits:
                    print("Half adder found for %s!" % filename)
                    print("Input island %d, output island %d, where" % (i+1, j+1))
                    print("angles %s correspond to %s." % (list(angles_i), input_bits))
                    print("#"*80)
                    found_halfadders.append((i, j, angles_i, input_bits))
    return found_halfadders


# TODO: write function that plots the low energies, given a certain input island.
# 
def plot_energy_levels(filename, input_island, trunc=5):
    '''
        Shows a matplotlib plot, with four lines in it. Each line corresponds to the
        input island being put at a certain input bit. For this input bit angle, a plot-line
        is then drawn corresponding to the energy of all possible configurations, in order.
        This will nicely illustrate whether the geometry is balanced, and how robust it is.
        You see, if one of these plot lines is too flat, it is easy for the geometry to go into
        this other, slightly higher energy state, which is of course undesirable.

        @param filename [str]: Relative path to mumax3 table.txt file.
        @param input_island [int]: The island that is considered the 'input island'.
        @param trunc [int] (5): Maximum amount of lowest-lying-energy-levels that are plotted/printed.
    '''
    mag_angles, energies, geom_angles = process_data(filename)
    y_trunc = 0
    for i in range(4):
        angle = geom_angles[input_island-1]+i*90
        collected_mag_angles, collected_energies = collect_orientation(input_island, angle, mag_angles, energies, geom_angles)
        if len(collected_energies) == 0:
            continue
        plt.plot(range(len(collected_energies)), collected_energies, label='Input %d°' % angle)
        y_trunc = max(y_trunc, max(collected_energies[:trunc+1]))
        print('Energy levels for input island %d at %d deg:' % (input_island, angle))
        print(collected_energies[:trunc+1])
        with np.printoptions(precision=0, suppress=True):
            print('Corresponding magnetization angles:')
            print(np.array(collected_mag_angles[:trunc+1]))
        print('-'*80)
    plt.title('Input island %d' % input_island)
    plt.legend()
    plt.xlabel('Energy level $n$')
    plt.ylabel('Energy [eV]')
    plt.grid(axis='x', color='grey', linestyle=':', linewidth=1)
    plt.xlim([0,trunc])
    plt.ylim([0,y_trunc])
    plt.gcf().tight_layout()
    plt.show()

if __name__ == "__main__":
    # TODO: Make the conversion auto-detect whether the folder is 'many_islands_interaction.out'.
    # convert_ovf('many_islands_interaction.out')
    # print('#'*80)

    check_if_halfadder('many_islands_interaction.out/table.txt')
    plot_energy_levels('many_islands_interaction.out/table.txt', 1)
    # get_lowest_energies('many_islands_interaction.out/table.txt', 3, verbose=False)
    # get_lowest_energies('many_islands_interaction.out/table.txt', 4, verbose=False)
    # get_lowest_energies('many_islands_interaction.out/table.txt', 6, verbose=False)
    # mag_angles, energies, geom_angles = process_data('many_islands_interaction.out/table.txt')
    # collect_orientation(1, 0, mag_angles, energies, geom_angles)

    # mag_angles, energies, geom_angles = process_data('attempts/table000006.txt')
    # collect_orientation(1, 0, mag_angles, energies, geom_angles)
    # get_lowest_energies('attempts/table000006.txt', 1, verbose=False)
    # check_if_halfadder('attempts/table000006.txt')
    # plot_energy_levels('attempts/table000006.txt', 1)
    # plot_energy_levels('attempts/table000010.txt', 3)
    # plot_energy_levels('attempts/table000012.txt', 1)

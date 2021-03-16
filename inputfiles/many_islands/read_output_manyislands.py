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
    for i in range(num_islands):
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


def collect_orientation(island, angle, mag_angles, energies, geom_angles, margin_mag=22.5, margin_energy=0.01, sort_energies=True):
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
        found_similar_entry = False
        for j, other in enumerate(mag_angles[:i]):
            if np.all(np.abs(anglediff(entry, other)) < margin_mag) and abs(energies[i]-energies[j]) < margin_energy: # True if the two entries are similar (both with angles and energy)
                found_similar_entry = True
                break
        if found_similar_entry:
            continue
        else:
            new_mag_angles.append(entry)
            new_energies.append(energies[i])
    return (new_mag_angles, new_energies)

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
            return (None, None)
        else:
            minEnergies.append(new_energies[0])
            minEnergyMagAngles.append(new_mag_angles[0])

    if not verbose:
        print('Lowest energies: \n', minEnergies)
        print('\nStable angles: \n', np.array(minEnergyMagAngles))
    return (np.array(minEnergyMagAngles), np.array(minEnergies))

def check_if_halfadder(filename):
    """
        For a given mumax3 table file at <filename>, this function looks if there is a half adder
        possible with the relaxed magnetizations from the table file.

        @param filename [str]: Relative path to mumax3 table.txt file.
        @returns: - False if no half adder was found.
                  - (i, j, angles_i, input_bits), where
                    i: Index of input island.
                    j: Index of output island.
                    angles_i: The four stable angles for the input island magnetization.
                    input_bits: The bit associated with each magnetization angle in angles_i.
    """
    _, _, geom_angles = process_data(filename)
    halfAdder = {0:0,1:1,2:1,3:2}
    for i in range(len(geom_angles)): # Input island <i>
        minEnergyMagAngles, _ = get_lowest_energies(filename, i+1)
        if minEnergyMagAngles is None: # Then no stable configurations exist for some angle of input <i> (e.g. because island <i> is fixed)
            continue
        angles_i = round_90(minEnergyMagAngles[:,i].T, delta=geom_angles[i])
        if np.all(np.array(angles_i) == angles_i[0]): # Then island <i> is fixed
            continue
        for j in range(len(geom_angles)): # Output island <j>
            if i == j or geom_angles[i] != geom_angles[j]: # If same island, or not same orientation, then skip
                continue
            angles_j = round_90(minEnergyMagAngles[:,j].T, delta=geom_angles[j])
            if np.all(np.array(angles_j) == angles_j[0]): # Then island <j> is fixed
                continue
            for perm in itertools.permutations([0,1,2,3]):
                d = {angles_i[n]:perm[n] for n in range(4)}
                input_bits = perm
                should_be_output_bits = list(map(halfAdder.get, input_bits))
                output_bits = list(map(d.get, angles_j))
                if should_be_output_bits == output_bits:
                    print("Half adder found!")
                    print("Input island %d, output island %d, where" % (i, j))
                    print("angles %s correspond to %s." % (list(angles_i), input_bits))
                    return (i, j, angles_i, input_bits)
    return False


# TODO: write function that plots the low energies, given a certain input island.
# This should be a matplotlib plot, with four lines in it. Each line corresponds to the
# input island being put at a certain input bit. For this input bit angle, a plot-line
# is then drawn corresponding to the energy of all possible configurations, in order.
# This will nicely illustrate whether the geometry is balanced, and how robust it is.
# You see, if one of these plot lines is too flat, it is easy for the geometry to go into
# this other, slightly higher energy state, which is of course undesirable.


if __name__ == "__main__":
    # TODO: Make the conversion auto-detect whether the folder is 'many_islands_interaction.out'.
    # convert_ovf('many_islands_interaction.out')
    # print('#'*80)

    # mag_angles, energies, geom_angles = process_data('many_islands_interaction.out/table.txt')
    # collect_orientation(1, 0, mag_angles, energies, geom_angles)
    # get_lowest_energies('many_islands_interaction.out/table.txt', 1)
    # check_if_halfadder('many_islands_interaction.out/table.txt')
    check_if_halfadder('attempts/table000006.txt')

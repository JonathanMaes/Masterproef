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
from itertools import product
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

def plot_probabilities(halfadderfile):
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
    collected_magangles_input = []
    collected_magangles_output = []
    largest_collection = 0 # This is the largest number of energy levels
    for i in range(4):
        new_mag_angles, new_energies = collect_orientation(input_island, geom_angles[input_island-1]+i*90, mag_angles, energies, geom_angles)
        collected_energies.append(new_energies)
        collected_magangles_input.append(new_mag_angles[:,input_island-1])
        collected_magangles_output.append(new_mag_angles[:,output_island-1])
        largest_collection = max(largest_collection, len(collected_energies))
    collected_magangles_input = [round_90(l, delta=geom_angles[input_island-1]) % 360 for l in collected_magangles_input]
    collected_magangles_output = [round_90(l, delta=geom_angles[output_island-1]) % 360 for l in collected_magangles_output]
    collected_inputs = [[halfadder_comb[int(angle//90)] for angle in l] for l in collected_magangles_input]
    collected_outputs = [[halfadder_comb[int(angle//90)] for angle in l] for l in collected_magangles_output]
    print(collected_energies)
    print(collected_inputs)
    print(collected_outputs)

    E_therm = 0.258 # [eV]

    portions = [np.exp(-line/E_therm) for line in collected_energies]
    portions_sum = [np.sum(line) for line in portions]
    # probabilities = [line/sum(portions_sum) for line in portions] # probability of (input,output) combination if all free to do what they want
    probabilities = [line/np.sum(line) for line in portions] # probability of output IF INPUT FIXED AT THE VALUE
    print(probabilities)

    flattened_energies = [i for line in collected_energies for i in line]
    flattened_inputs = [i for line in collected_inputs for i in line]
    flattened_outputs = [i for line in collected_outputs for i in line]
    flattened_inout = [(w, flattened_outputs[i]) for i, w in enumerate(flattened_inputs)]
    # print(flattened_inout)

    correct = [(0,0), (1,1), (2,1), (3,2)]

    # plot as function of thermal energy
    E_therms = np.arange(0, 1, 0.001)
    portionss = np.array([[np.exp(-E/kBT) for kBT in E_therms] for E in flattened_energies])
    probabilitiess = [[w/np.sum(portionss[:,i]) for i,w in enumerate(comb)] for comb in portionss]
    for i, line in enumerate(probabilitiess):
        plt.plot(E_therms, line, color='b' if flattened_inout[i] in correct else 'r')
    plt.xlabel('Thermal energy [eV]')
    plt.ylabel('Probability [eV]')
    plt.show()

    flattened_inputfixed_probs = [i for line in [line/np.sum(line) for line in portions] for i in line]
    print(flattened_inputfixed_probs)

    # Let's now propagate through two half adders
    all_states = []
    all_probs = []
    for A in product([0,1], repeat=3):
        A1, A2, A3 = A

        possible_states = [comb for comb in product([0,1], repeat=4)] # (D, S3, S2, S1)
        for state in possible_states:
            D, S3, S2, S1 = state

            prob = 1

            # What is the chance that the first half adder spits out (S1, D)?
            out1 = S1 + 2*D
            in1 = A1 + 2*A2
            if (in1, out1) in flattened_inout:
                idx1 = flattened_inout.index((in1, out1))
                prob *= flattened_inputfixed_probs[idx1]
            else:
                continue # Combination is impossible
                
            # What is the chance that second half adder then spits out (S2, S3)?
            out2 = S2 + 2*S3
            in2 = D + 2*A3
            if (in2, out2) in flattened_inout:
                idx2 = flattened_inout.index((in2, out2))
                prob *= flattened_inputfixed_probs[idx2]
            else:
                continue # Combination is impossible
        
            all_states.append((*A, *state))
            all_probs.append(prob)
    print(all_probs)
    
    # Now we have to collect by output, eliminating the intermediate D
    sifted_states = []
    sifted_probs = []
    for i, thing in enumerate(all_states):
        sifted_state = (*thing[0:3], *thing[4:])
        if sifted_state not in sifted_states:
            sifted_states.append(sifted_state)
            sifted_probs.append(all_probs[i])
        else:
            idx = sifted_states.index(sifted_state)
            sifted_probs[idx] += all_probs[i]
    print(sifted_states)
    print(sifted_probs)
    print(max(sifted_probs))

    image_array = np.zeros((2**3, 2**3))
    input_combs = [i for i in product([0,1], repeat=3)]
    output_combs = [i for i in product([0,1], repeat=3)]
    for i, input_comb in enumerate(input_combs):
        for j, output_comb in enumerate(output_combs):
            state = (*input_comb, *output_comb)
            if state in sifted_states:
                idx = sifted_states.index(state)
                prob = sifted_probs[idx]
                image_array[i, j] = prob
    
    
    fig = plt.figure(figsize=(10.0, 5.0))
    ax = fig.add_subplot(111)
    im = ax.imshow(image_array, cmap=cm.get_cmap('Greys'))
    plt.colorbar(im)
    plt.scatter([0, 2, 1, 3, 1, 3, 2, 4], [0, 1, 2, 3, 4, 5, 6, 7])
    ax.set_xticks(np.arange(image_array.shape[1]))
    ax.set_yticks(np.arange(image_array.shape[0]))
    ax.set_xlabel('Output')
    ax.set_ylabel('Input')
    ax.set_xticklabels(output_combs)
    ax.set_yticklabels(input_combs)
    ax.tick_params(top=True, bottom=False,
                   labeltop=True, labelbottom=False)

    # Rotate the tick labels and set their alignment.
    plt.setp(ax.get_xticklabels(), rotation=-30, ha="right",
             rotation_mode="anchor")

    plt.gcf().tight_layout()

    plt.show()


        






    # # plt.gcf().subplots_adjust(bottom=0.15, left=0.15, right=0.95)
    # plt.gcf().tight_layout()
    # plt.savefig(os.path.join(outDir, os.path.split(halfadderfile)[1].replace('.txt', '_energylevels.pdf')))
    # plt.show()




if __name__ == "__main__":
    pass
    plot_probabilities('Results/Halfadder000006/table(d170,s-60).txt')
    # plot_probabilities('Results/Halfadder000006/tableside(d50,s130).txt')

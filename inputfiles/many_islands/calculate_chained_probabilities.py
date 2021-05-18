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

def plot_probabilities(halfadderfile, logicconfig=1):
    """
        @param halfadderfile [str]: The mumax table.txt containing all possible states of the half adder.
        @param logicconfig [int] (1): The configuration in which half adders are chained together.
            1: Two half adders in series, with the second half adder taking the carry of the first.
            2: Add two two-bit numbers.
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
    print('Energies:', collected_energies)
    print('Inputs:  ', collected_inputs)
    print('Outputs: ', collected_outputs)

    E_therm = 0.0258 # [eV]

    flattened_energies = [i for line in collected_energies for i in line]
    flattened_inputs = [i for line in collected_inputs for i in line]
    flattened_outputs = [i for line in collected_outputs for i in line]
    flattened_inout = [(w, flattened_outputs[i]) for i, w in enumerate(flattened_inputs)]


    # # plot as function of thermal energy
    # E_therms = np.arange(0, 1, 0.001)
    # portionss = np.array([[np.exp(-E/kBT) for kBT in E_therms] for E in flattened_energies])
    # probabilitiess = [[w/np.sum(portionss[:,i]) for i,w in enumerate(comb)] for comb in portionss]
    # correct = [(0,0), (1,1), (2,1), (3,2)]
    # for i, line in enumerate(probabilitiess):
    #     plt.plot(E_therms, line, color='b' if flattened_inout[i] in correct else 'r')
    # plt.xlabel('Thermal energy [eV]')
    # plt.ylabel('Probability [eV]')
    # plt.show()



    # Let's now propagate through two half adders
    if logicconfig == 1:
        fig = plt.figure(figsize=(6.5, 5.0))
        ax = fig.add_subplot(111)

        all_states = []
        all_energies = []
        all_probs = []
        for A in product([0,1], repeat=3):
            A1, A2, A3 = A
            current_energies = []

            possible_states = [comb for comb in product([0,1], repeat=4)] # (D, S1, S2, S3)
            for state in possible_states:
                D, S1, S2, S3 = state

                energy = 0

                # What is the chance that the first half adder spits out (S1, D)?
                out1 = S1 + 2*D
                in1 = A1 + 2*A2
                if (in1, out1) in flattened_inout:
                    idx1 = flattened_inout.index((in1, out1))
                    energy += flattened_energies[idx1]
                else:
                    continue # Combination is impossible
                    
                # What is the chance that second half adder then spits out (S2, S3)?
                out2 = S2 + 2*S3
                in2 = D + 2*A3
                if (in2, out2) in flattened_inout:
                    idx2 = flattened_inout.index((in2, out2))
                    energy += flattened_energies[idx2]
                else:
                    continue # Combination is impossible
            
                all_states.append((*A, *state))
                all_energies.append(energy)
                current_energies.append(energy)
        
            relative_probs = [np.exp(-E/E_therm) for E in current_energies]
            all_probs += [i/sum(relative_probs) for i in relative_probs]


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
        ax.set_xticklabels(output_combs)
        ax.set_yticklabels(input_combs)
                        
        plt.scatter([0, 2, 4, 6, 4, 6, 2, 1], [0, 1, 2, 3, 4, 5, 6, 7])
    
    
    if logicconfig == 2:
        fig = plt.figure(figsize=(4.5, 5.0))
        ax = fig.add_subplot(111)
        
        all_states = []
        all_energies = []
        all_probs = []
        for A in product([0,1], repeat=4):
            A2, A1, B2, B1 = A
            current_energies = []

            possible_states = [comb for comb in product([0,1], repeat=7)] # (C1, S0, C2, C3, S3, S2, S1)
            for state in possible_states:
                C1, S0, C2, C3, S3, S2, S1 = state

                energy = 0

                # What is the chance that the first half adder spits out (S1, C1)?
                out1 = S1 + 2*C1
                in1 = A1 + 2*B1
                if (in1, out1) in flattened_inout:
                    idx1 = flattened_inout.index((in1, out1))
                    energy += flattened_energies[idx1]
                else:
                    continue # Combination is impossible
                    
                # What is the chance that third half adder spits out (S0, C3)?
                out2 = S0 + 2*C3
                in2 = A2 + 2*B2
                if (in2, out2) in flattened_inout:
                    idx2 = flattened_inout.index((in2, out2))
                    energy += flattened_energies[idx2]
                else:
                    continue # Combination is impossible
                
                # What is the chance that the second half adder then spits out (S2, C2)?
                out1 = S2 + 2*C2
                in1 = C1 + 2*S0
                if (in1, out1) in flattened_inout:
                    idx1 = flattened_inout.index((in1, out1))
                    energy += flattened_energies[idx1]
                else:
                    continue # Combination is impossible

                # Assume ideal or gate spitting out S3
                if (C2 or C3) != S3:
                    continue # Combination is impossible
            
                all_states.append((*A, *state))
                all_energies.append(energy)
                current_energies.append(energy)
        
            relative_probs = [np.exp(-E/E_therm) for E in current_energies]
            all_probs += [i/sum(relative_probs) for i in relative_probs]
        
        # Now we have to collect by output, eliminating the intermediate D
        sifted_states = []
        sifted_probs = []
        for i, thing in enumerate(all_states):
            sifted_state = (*thing[0:4], *thing[8:])
            if sifted_state not in sifted_states:
                sifted_states.append(sifted_state)
                sifted_probs.append(all_probs[i])
            else:
                idx = sifted_states.index(sifted_state)
                sifted_probs[idx] += all_probs[i]

        image_array = np.zeros((2**4, 2**3))
        input_combs = [i for i in product([0,1], repeat=4)]
        output_combs = [i for i in product([0,1], repeat=3)]
        for i, input_comb in enumerate(input_combs):
            for j, output_comb in enumerate(output_combs):
                state = (*input_comb, *output_comb)
                if state in sifted_states:
                    idx = sifted_states.index(state)
                    prob = sifted_probs[idx]
                    image_array[i, j] = prob
        ax.set_xticklabels([4*S3+2*S2+1*S1 for S3, S2, S1 in output_combs])
        ax.set_yticklabels(['%d+%d' % (2*A2+A1, 2*B2+B1) for A2, A1, B2, B1 in input_combs])

        plt.scatter([0, 1, 2, 3, 1, 2, 3, 4, 2, 3, 4, 5, 3, 4, 5, 6], np.arange(16), s=9)
    
    
    im = ax.imshow(image_array, vmin=0, cmap=cm.get_cmap('Greys'))
    cbar = plt.colorbar(im)
    cbar.set_label(r'Forward probability', rotation=270, labelpad=25)
    ax.set_xticks(np.arange(image_array.shape[1]))
    ax.set_yticks(np.arange(image_array.shape[0]))
    ax.set_xlabel('Output')
    ax.xaxis.set_label_position('top') 
    ax.set_ylabel('Input')
    ax.tick_params(top=True, bottom=False,
                   labeltop=True, labelbottom=False)

    # Rotate the tick labels and set their alignment, if necessary.
    if logicconfig in [1]:
        plt.setp(ax.get_xticklabels(), rotation=90, rotation_mode='anchor', horizontalalignment='left', verticalalignment='center')

    plt.gcf().tight_layout()
    plt.savefig(halfadderfile.replace('.txt', '_forward%d_probabilities_T%seV.pdf' % (logicconfig, E_therm)))

    plt.show()


if __name__ == "__main__":
    pass
    plot_probabilities('Results/Halfadder000006/table(d170,s-60).txt', logicconfig=1)
    # plot_probabilities('Results/Halfadder000006/table(d170,s-60).txt', logicconfig=2)
    # plot_probabilities('Results/Halfadder000006/tableside(d50,s130).txt', logicconfig=1)
    # plot_probabilities('Results/Halfadder000006/tableside(d50,s130).txt', logicconfig=2)

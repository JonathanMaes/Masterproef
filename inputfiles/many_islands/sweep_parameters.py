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
font = {'size':16}
matplotlib.rc('font', **font)

from read_output_manyislands import check_if_halfadder, read_mumax3_table
from generate_mumax3_inputfile import generate_mumax3_inputfile, Island

# Define some filenames
attemptcode = 'test'
# attemptcode = '000006'
outDir = 'Results/Sweeps/Sweep_%s' % attemptcode
if not os.path.exists(outDir):
    os.makedirs(outDir)

def sweep():
    # Because it is too complicated to make this function accept arguments, just change things manually
    ########* CHANGE INPUT VARIABLES BELOW *########
    var1, var1_max, var1_step, var1_name, var1_unit = -160, -200, -10, 'd', 'nm'
    var2, var2_max, var2_step, var2_name, var2_unit = 0, 50, 10, 's', 'nm'
    #! WHEN CHANGING THE MEANING OF var1 OR var2, DONT FORGET TO ADJUST THE generate_mumax3_inputfile CALL IN THE LOOP ACCORDINGLY
    #! Also dont forget to make sure the islands are still on the right grid positions when sweeping these variables.
    ########* CHANGE INPUT VARIABLES ABOVE *########
    outTableTxt = os.path.join(outDir, 'table(%s,%s).txt' % (var1_name, var2_name))

    

    # Define some ranges and iterate over var1 and var2
    var1_range = list(np.arange(var1, var1_max+var1_step/2, var1_step))
    var2_range = list(np.arange(var2, var2_max+var2_step/2, var2_step))
    for i, var1 in enumerate(var1_range):
        for j, var2 in enumerate(var2_range):
            print('Testing %s=%s%s, %s=%s%s...' % (var1_name, var1, var1_unit, var2_name, var2, var2_unit))
            # Create mumax3 inputfile to test half adder for these values of var1, var2
            generate_mumax3_inputfile(2, [ # Basic half adder of type 000006
                    Island(-128, 0, 0),
                    Island(0, 0, 0),
                    Island(var2, var1, pi/2, fixed=True)
                ], rho=0.66, L=100, extra_columns={'var1-%s [%s]' % (var1_name, var1_unit): var1, 'var2-%s [%s]' % (var2_name, var2_unit): var2})
            # Run the newly generated mumax3 inputfile
            os.system('mumax3 many_islands_interaction.mx3')
            # Save this table.txt and add it to a huge table.txt in the outDir folder
            if i == j == 0:
                shutil.copyfile('many_islands_interaction.out/table.txt', outTableTxt)
            else:
                with open('many_islands_interaction.out/table.txt', 'r') as tabletxt:
                    text = [line for line in tabletxt]
                with open(outTableTxt, 'a') as outTable:
                    for line in text[1:]: # First line is the header, which we dont want to copy
                        outTable.write(line)


def plot_sweep(sweepfile):
    """
        @param sweepfile [str]: The relative path to the "table(var1,var2).txt".
    """
    table = read_mumax3_table(sweepfile)
    with open(sweepfile) as tabletxt:
        headers = tabletxt.readline().split('\t')
        for header in headers:
            if header.startswith('var1-'):
                var1_name, var1_unit = header.split('var1-')[1].split(' ')
                var1_unit = var1_unit.strip()[1:-1] # remove parentheses at ends of string
            elif header.startswith('var2-'):
                var2_name, var2_unit = header.split('var2-')[1].split(' ')
                var2_unit = var2_unit.strip()[1:-1] # remove parentheses at ends of string

    var1_range = [val for val, _ in table.groupby("var1-%s" % var1_name, sort=True)]
    var2_range = [val for val, _ in table.groupby("var2-%s" % var2_name, sort=True)]
    var1_step = var1_range[1] - var1_range[0] if len(var1_range) > 1 else var1_range[0]
    var2_step = var2_range[1] - var2_range[0] if len(var2_range) > 1 else var2_range[0]
    is_halfadder_grid = np.zeros((len(var1_range), len(var2_range))) # List of lists, main list enumerates var1, sublists contain var2 enumeration
    for i, subset in enumerate(table.groupby("var1-%s" % var1_name, sort=True)):
        subtable = subset[1] # Element 0 is the value of var1 in this loop
        for j, subsubset in enumerate(subtable.groupby("var2-%s" % var2_name, sort=True)):
            subsubtable = subsubset[1] # Element 0 is the value of var2 in this loop
            halfadders = check_if_halfadder(subsubtable, verbose=False)
            if len(halfadders) != 0:
                print(halfadders[0][3])
                is_halfadder_grid[i][j] = 1

    # Plot result
    fig = plt.figure(figsize=(7.0, 5.0))
    lim_1 = np.array([var1_range[0]-var1_step/2, var1_range[-1]+var1_step/2])
    lim_2 = np.array([var2_range[0]-var2_step/2, var2_range[-1]+var2_step/2])

    ax = fig.add_subplot(111)
    im = ax.imshow(np.matrix(is_halfadder_grid)[:,::-1].T, vmin=0, vmax=1, interpolation='nearest', cmap=cm.get_cmap('inferno'), extent=[lim_1[0], lim_1[1], lim_2[0], lim_2[1]])
    ax.set_aspect('auto') # Stretch figure to fit pdf nicely.
    plt.xlabel(var1_name + ' [%s]' % var1_unit)
    plt.ylabel(var2_name + ' [%s]' % var2_unit)
    cbar = fig.colorbar(im)
    cbar.set_label('Is half adder?', rotation=270, labelpad=25)

    # plt.xticks(var1_range)
    # plt.yticks(var2_range)

    plt.gcf().subplots_adjust(bottom=0.15, left=0.15, right=0.95)
    plt.gcf().tight_layout()
    plt.savefig(os.path.join(outDir, os.path.split(sweepfile)[1].replace('.txt', '.pdf')))
    plt.show()


if __name__ == "__main__":
    plot_sweep('Results/Sweeps/Sweep_test/table(d130-220_10,s-100-100_10).txt')
    # plot_sweep('Results/Sweeps/Sweep_test/table(d100-210_2,s20).txt')

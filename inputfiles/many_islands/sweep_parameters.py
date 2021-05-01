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
from matplotlib.colors import ListedColormap
font = {'size':16}
matplotlib.rc('font', **font)

from read_output_manyislands import check_if_halfadder, read_mumax3_table, process_data, collect_orientation
from generate_mumax3_inputfile import generate_mumax3_inputfile, Island

# Define some filenames
# attemptcode = 'test'
attemptcode = '000006'
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

    
def contour_rect(im, extent):
    """ This code is not very efficient, but since the imshows are quite small this is not a problem.
        To plot these lines, use

        for line in lines:
            plt.plot(line[1], line[0], **kwargs)

        @param im [list of lists]: The grid that is plotted in the imshow.
        @param extent [list(4)]: List containing [xmin, xmax, ymin, ymax]. This is the same as the imshow 'extent' kwarg.
        @return [list]: List of tuples of the form ([y0, y1], [x0, x1]), representing lines from (x0, y0) to (x1, y1).
    """
    xmin = extent[0]
    xmax = extent[1]
    ymin = extent[2]
    ymax = extent[3]
    ny, nx = im.shape
    dy, dx = (ymax - ymin)/ny, (xmax - xmin)/nx

    pad = np.pad(im, [(1, 1), (1, 1)], mode='constant')  # zero padding
    im0 = np.abs(np.diff(pad, n=1, axis=0))[:, 1:] # Diff in y-direction
    im1 = np.abs(np.diff(pad, n=1, axis=1))[1:, :] # Diff in x-direction

    lines = []
    for ii, jj in np.ndindex(im0.shape):
        # Format: lines += [([y0, y1], [x0, x1])] for line from (x0, y0) to (x1, y1)
        if im0[ii, jj] > 0: 
            lines += [([ymin + dy*ii, ymin + dy*ii], [xmin + dx*jj, xmin + dx*(jj+1)])]
        if im1[ii, jj] > 0:
            lines += [([ymin + dy*ii, ymin + dy*(ii+1)], [xmin + dx*jj, xmin + dx*jj])]
    return lines

def replace_with_dict(ar, dic):
    k = np.array(list(dic.keys()))
    v = np.array(list(dic.values()))
    sidx = k.argsort()
    return v[sidx[np.searchsorted(k,ar,sorter=sidx)]]

def plot_sweep(sweepfile, swap_axes=False, do=('types', 'balanced1')):
    """
        @param sweepfile [str]: The relative path to the "table(var1,var2).txt".
        @param swap_axes [bool] (False): If True, var1 is plotted on the y-axis and var2 on the x-axis.
                                         If False, var1 on the x-axis and var2 on the y-axis, as normal.
        @param do [tuple] ('types', 'balanced1'): All the actions in this tuple are plotted and shown. These actions are:
            - 'types': Imshow where the image has several colors, each color corresponding to one type of halfadder.
                       A type of halfadder means which angles correspond to which quaternary number, e.g. (0, 2, 1, 3)
            - 'balanced1': Imshow where min_a(E_{a,1}) - max_a(E_{a,0}) is plotted (diff highest ground state with lowest first excited state).
            - 'balanced2': Imshow where max_a(E_{a,0}) - min_a(E_{a,0}) is plotted (diff highest ground state with lowest ground state). 
    """
    # Things this function could do:
    # - Plot the different half adder types with colors and colorbar, to know what each contour means in the other types of plots
    # - Plot the balancedness type 1: min_a(E_{a,1}) - max_a(E_{a,0})
    #                         type 2: max_a(E_{a,0}) - min_a(E_{a,0})

    if swap_axes:
        var1_colstart = 'var2-'
        var2_colstart = 'var1-'
    else:
        var1_colstart = 'var1-'
        var2_colstart = 'var2-'
    
    #### Read table and determine all relevant half adder quantities
    table = read_mumax3_table(sweepfile)
    with open(sweepfile) as tabletxt:
        headers = tabletxt.readline().split('\t')
        for header in headers:
            if header.startswith(var1_colstart):
                var1_name, var1_unit = header.split(var1_colstart)[1].split(' ')
                var1_unit = var1_unit.strip()[1:-1] # remove parentheses at ends of string
            elif header.startswith(var2_colstart):
                var2_name, var2_unit = header.split(var2_colstart)[1].split(' ')
                var2_unit = var2_unit.strip()[1:-1] # remove parentheses at ends of string

    var1_range = [val for val, _ in table.groupby(var1_colstart + var1_name, sort=True)]
    var2_range = [val for val, _ in table.groupby(var2_colstart + var2_name, sort=True)]
    var1_sweeped, var2_sweeped = len(var1_range) > 1, len(var2_range) > 1
    var1_step = var1_range[1] - var1_range[0] if var1_sweeped else var2_range[1] - var2_range[0]
    var2_step = var2_range[1] - var2_range[0] if var2_sweeped else var1_range[1] - var1_range[0]

    is_halfadder_grid = np.zeros((len(var1_range), len(var2_range))) # List of lists, main list enumerates var1, sublists contain var2 enumeration
    halfadder_types = []
    halfadder_type_grid = np.zeros_like(is_halfadder_grid) # Zero if no half adder, 1 corresponds to halfadder_types[0], 2 to halfadder_types[1], ...
    balance1_grid = np.zeros_like(is_halfadder_grid)
    balance1_grid[:] = np.nan
    balance2_grid = np.zeros_like(is_halfadder_grid)
    balance2_grid[:] = np.nan
    for i, subset in enumerate(table.groupby(var1_colstart + var1_name, sort=True)):
        subtable = subset[1] # Element 0 is the value of var1 in this loop
        for j, subsubset in enumerate(subtable.groupby(var2_colstart + var2_name, sort=True)):
            subsubtable = subsubset[1] # Element 0 is the value of var2 in this loop
            halfadders = check_if_halfadder(subsubtable, verbose=False)
            print(halfadders)
            if len(halfadders) != 0:
                # Do everything that has to do with being a half adder or not, and the type of halfadder (which is input? angles in which order?)
                input_idx, output_idx, angles_i, input_bits = halfadders[0]
                input_bits = [x for _, x in sorted(zip(angles_i, input_bits))] # Sort input bits in the same manner as angles_i
                angles_i = sorted(angles_i)
                input_bits = input_bits[2:] + input_bits[:2] # Make angles start at or above 0 instead of -180°
                angles_i = [angle + 180 for angle in angles_i] # Make angles start at or above 0 instead of -180°
                is_halfadder_grid[i][j] = 1
                halfadder_type = (input_idx, output_idx, *input_bits)
                if halfadder_type not in halfadder_types: # Element 3 contains which angles correspond to which number e.g. (0, 2, 1, 3)
                    halfadder_types.append(halfadder_type)
                idx = halfadder_types.index(halfadder_type)
                halfadder_type_grid[i][j] = idx+1

                # Do everything that has to do with energies
                mag_angles, energies, geom_angles = process_data(subsubtable)
                input_island = input_idx
                all_energies = []
                largest_collection = 0
                for a in range(4):
                    angle = geom_angles[input_island-1]+a*90
                    collected_mag_angles, collected_energies = collect_orientation(input_island, angle, mag_angles, energies, geom_angles)
                    all_energies.append(collected_energies)
                    largest_collection = max(largest_collection, len(collected_energies))
                all_energies = np.array([list(l) + [np.infty]*(largest_collection - len(l)) for l in all_energies])
                balance1_grid[i][j] = np.min(all_energies[:,1]) - np.max(all_energies[:,0])
                balance2_grid[i][j] = np.max(all_energies[:,0]) - np.min(all_energies[:,0])
                # balance1_grid[i][j] = np.max()
    
    n_types = len(halfadder_types)
    # Sort the types numbering so they are not in a weird order. This only affects the behind-the-scenes numbering of the types, not the halfadders themselves.
    sorted_types = sorted(halfadder_types)
    sorted_types_dict = {0:0, **{i+1:sorted_types.index(typ)+1 for i, typ in enumerate(halfadder_types)}}
    halfadder_type_grid = replace_with_dict(halfadder_type_grid, sorted_types_dict)
    halfadder_types = sorted_types

    #### Plots
    # First, prepare the halfadder type contours (to make the contours nicely follow the pixel shape, the image is first upscaled by a factor 64)
    lim_x = np.array([var1_range[0]-var1_step/2, var1_range[-1]+var1_step/2])
    lim_y = np.array([var2_range[0]-var2_step/2, var2_range[-1]+var2_step/2])
    extent = [lim_x[0], lim_x[1], lim_y[0], lim_y[1]]
    contourlines = contour_rect(halfadder_type_grid.transpose(), extent)
    print('Half adder types in order (input, output, %d deg, %d deg, %d deg, %d deg): %s' % (*angles_i, halfadder_types))
    def draw_contour(ax, **kwargs):
        for line in contourlines:
            ax.plot(line[1], line[0], **kwargs)
    

    ## Plot the type of half adder e.g. (0, 2, 1, 3)
    if 'types' in do:
        fig = plt.figure(figsize=(7.0, 5.0))
        ax = fig.add_subplot(111)

        plot_grid = np.transpose(halfadder_type_grid)
        cMapDiscr = ListedColormap(['white', '#5555ff', 'red', 'orange', 'darkgreen'][0:n_types+1])
        im = ax.imshow(plot_grid, vmin=0, vmax=n_types+1, origin='lower', interpolation='nearest', cmap=cMapDiscr, extent=extent) #, vmin=0, vmax=1
        cbar = fig.colorbar(im)
        # cbar.set_label('Is half adder?', rotation=270, labelpad=25)
        cbar.ax.get_yaxis().set_ticks([0.5+i for i in range(n_types+1)])
        cbar.ax.get_yaxis().set_ticklabels(['No half adder'] + ['In %d %s' % (tup[0], tup[2:]) for tup in halfadder_types])
        # for j, lab in enumerate(halfadder_types): # Put half adder type INSIDE the colorbar, but then text is quite smol
        #     cbar.ax.text(1.5, 1.5 + j, lab, ha='center', va='center', rotation=270, size=10)

        if var1_sweeped and var2_sweeped: # Stretch figure to fit pdf nicely if both vars sweeped, but keep square pixels otherwise.
            ax.set_aspect('auto')
        if var1_sweeped:
            plt.xlabel(var1_name + ' [%s]' % var1_unit)
        else:
            plt.xticks([])
        if var2_sweeped:
            plt.ylabel(var2_name + ' [%s]' % var2_unit)
        else:
            plt.yticks([])
        
        draw_contour(ax, color='black', alpha=1)

        plt.gcf().subplots_adjust(bottom=0.15, left=0.15, right=0.95)
        plt.gcf().tight_layout()
        plt.savefig(os.path.join(outDir, os.path.split(sweepfile)[1].replace('.txt', '_types.pdf')))
        plt.show()
    
    ## Plot the balancedness
    for i in [1,2]:
        if f'balanced{i}' in do:
            fig = plt.figure(figsize=(7.0, 5.0))
            ax = fig.add_subplot(111)

            plot_grid = np.transpose(eval(f"balance{i}_grid")) # eval is bad practice but idc
            im = ax.imshow(plot_grid,  origin='lower', interpolation='nearest', cmap=cm.get_cmap('inferno'), extent=extent) #, vmin=0, vmax=n_types+1,
            cbar = fig.colorbar(im)
            cbar.set_label(f'Balanced? (rule {i})', rotation=270, labelpad=25)

            if var1_sweeped and var2_sweeped: # Stretch figure to fit pdf nicely if both vars sweeped, but keep square pixels otherwise.
                ax.set_aspect('auto')
            if var1_sweeped:
                plt.xlabel(var1_name + ' [%s]' % var1_unit)
            else:
                plt.xticks([])
            if var2_sweeped:
                plt.ylabel(var2_name + ' [%s]' % var2_unit)
            else:
                plt.yticks([])
            
            draw_contour(ax, color='black', alpha=1)

            plt.gcf().subplots_adjust(bottom=0.15, left=0.15, right=0.95)
            plt.gcf().tight_layout()
            plt.savefig(os.path.join(outDir, os.path.split(sweepfile)[1].replace('.txt', f'_balanced{i}.pdf')))
            plt.show()



if __name__ == "__main__":
    # plot_sweep('Results/Sweeps/Sweep_000006/table(d100-220_10,s-100-100_10).txt', swap_axes=True, do=('types'))
    # plot_sweep('Results/Sweeps/Sweep_000006/table(d100-220_10,s-100-100_10).txt', swap_axes=True, do=('balanced1'))
    plot_sweep('Results/Sweeps/Sweep_000006/table(d100-220_10,s-100-100_10).txt', swap_axes=True, do=('balanced2'))
    # plot_sweep('Results/Sweeps/Sweep_000006/table(d100-210_2,s20).txt')

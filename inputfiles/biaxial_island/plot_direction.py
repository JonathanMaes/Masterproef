"""
    PLOTS: - 'Relaxed magnetization angle' AS FUNCTION OF 'Time'
    Variable 'FANCY': if True, angles are plotted in a continuous manner past [-180°, 180°].
                      if False, angles are constrained to [-180°, 180°].
    Variable 'groupBy': can be used to group by occurences of one table column (doesn't work entirely).
"""
import math
import matplotlib
import os
import numpy as np
import matplotlib.pyplot as plt
import pandas as pd

font = {'size':16}
matplotlib.rc('font', **font)


def read_mumax3_table(filename):
    """Puts the mumax3 output table in a pandas dataframe"""
    
    table = pd.read_csv(filename, sep='\t')
    table.columns = ' '.join(table.columns).split()[1::2]
    
    return table

def plot(filename, save=False, inset=None, groupBy=None):
    '''
        @param filename [str]: relative path of the mumax table
        @param save [bool] <False>: Whether or not to save the created plot.
        @param inset [list(4 or 6)]: x0, y0, dx, dy(, Delta_t, dt) of inset. If not specified, no inset is generated.
            x0, y0, dx, dy specify the physical location and size of the inset on the larger plot, as a ratio.
            Delta_t [ns] (optional)  is the center time value of the inset plot. If not specified, this is the middle of the main plot.
            dt [ns] (optional) is half of the time range of the inset plot. If not specified, this is 5ns.
        @param groupBy [str]: If specified, a plot will be made for each separate value in the column named <groupBy>.
    '''
    table = read_mumax3_table(filename)

    FANCY = True # If FANCY: follow angles past [-180°, 180°]
    INSET = bool(inset)
    if groupBy: # If groupBy is a string: multiple plots grouped by that property
        subsets = [subset[1] for subset in table.groupby(groupBy)]
        legend = ['%s: %s' % (groupBy, subset[0]) for subset in table.groupby(groupBy)]
    else:
        subsets = [table]

    _, ax = plt.subplots(figsize=[8.0, 5.0]) # First element is 'fig' but is not used
    if INSET:
        axins = ax.inset_axes(inset[:4])
    angleRange = [0, 0]
    timeRange = [0, 0]
    for subtable in subsets:
        angles = np.arctan2(subtable["my"], subtable["mx"])*180/math.pi
        timeRange[0] = min(timeRange[0], np.min(subtable["t"])*1e9)
        timeRange[1] = max(timeRange[1], np.max(subtable["t"])*1e9)
        if not FANCY:
            ax.plot(subtable["t"]*1e9, angles)
            if INSET:
                axins.plot(subtable["t"]*1e9, angles)
        else:
            previousAngles = np.array(angles)[:-1]
            nextAngles = np.array(angles)[1:]
            offsets = (np.abs(previousAngles - nextAngles) > 180)*((previousAngles > nextAngles)*2 - 1)*360
            offset = np.cumsum(offsets)
            fancyAngles = np.append([angles[0]], angles[1:] + offset)
            angleRange[0] = min(angleRange[0], np.min(fancyAngles))
            angleRange[1] = max(angleRange[1], np.max(fancyAngles))
            ax.plot(subtable["t"]*1e9, fancyAngles)
            if INSET:
                axins.plot(subtable["t"]*1e9, fancyAngles)

    if groupBy:
        plt.legend(legend)

    plt.grid(axis='y', color='grey', linestyle=':', linewidth=1)
    plt.xlabel(r'$t$ [ns]')
    plt.ylabel(r'angle [°]')
    plt.xlim(timeRange[0], timeRange[1])
    angle_interval = 90
    yticks = np.arange(angleRange[0]//angle_interval*angle_interval, angleRange[1]//angle_interval*angle_interval+angle_interval+1, angle_interval)
    plt.yticks(yticks)
    if INSET:
        n = len(subtable["t"])
        t_step = (subtable["t"].iloc[-1] - subtable["t"][0])*1e9/n
        dt, Delta_n = 5, n//2
        Delta_t = Delta_n*t_step
        if len(inset) > 4:
            Delta_t = inset[4]
            Delta_n = int(Delta_t//t_step)
            if len(inset) > 5:
                dt = inset[5]
        if min(subtable["t"]*1e9) > (Delta_t-dt) or max(subtable["t"]*1e9) < (Delta_t + dt):
            print(dt, Delta_t)
            raise ValueError('Inset plot time range not fully in time range of inputfile.')
        dn = int(dt//t_step)
        x1, x2 = subtable["t"][Delta_n-dn]*1e9, subtable["t"][Delta_n+dn]*1e9
        y1, y2 = min(fancyAngles[Delta_n-dn:Delta_n+dn])-10, max(fancyAngles[Delta_n-dn:Delta_n+dn])+10
        axins.set_xlim(x1, x2)
        axins.set_ylim(y1, y2)
        axins.set_xticks([x1, subtable["t"][Delta_n]*1e9, x2])
        axins.set_xticklabels([int(round(x1)), int(round(subtable["t"][Delta_n]*1e9)), int(round(x2))], fontsize=10)
        axins.set_yticks(yticks[(y1 < yticks)*(yticks < y2)])
        axins.set_yticklabels([int(round(i)) for i in yticks[(y1 < yticks)*(yticks < y2)]], fontsize=10)
        axins.grid(axis='y', color='grey', linestyle=':', linewidth=1)
        ax.indicate_inset_zoom(axins, zorder=10, alpha=0.8)
    
    plt.gcf().tight_layout()
    
    if save:
        outDir = 'Figures/Switching'
        if not os.path.exists(outDir):
            os.makedirs(outDir)
        outFileName = os.path.join(outDir, os.path.splitext(os.path.basename(filename).split('table_')[-1])[0]) + ('_INSET.pdf' if INSET else '.pdf')
        plt.savefig(outFileName)

    plt.show()


if __name__ == "__main__":
    # plot('biaxial_island_switching_plus.out/table_65x100_300K_alpha0.01_1µs_4nm.txt', save=1, inset=[0.6, 0.6, 0.37, 0.37, 470, 5])
    plot('biaxial_island_switching_plus.out/table_65x100_300K_alpha0.1_1µs_4nm.txt', save=1, inset=[0.6, 0.1, 0.37, 0.37, 420, 5])
    # inFileName = 'biaxial_island_switching_plus.out/table_65x100_300K_alpha0.1_1µs_4nm.txt'
    # inFileName = 'biaxial_island_switching_plus.out/table_65x100_300K_alpha0.01_1µs_4nm.txt'
    # inFileName = 'biaxial_island_switching_plus.out/table_65x100_350K_alpha0.01_1µs_4nm.txt'
    # inFileName = 'biaxial_island_switching_plus.out/table_65x100_273K_alpha0.01_1µs_4nm.txt'
    # inFileName = 'biaxial_island_switching_plus.out/table_100x100_350K_alpha0.01_1µs_4nm.txt'
    # inFileName = 'biaxial_island_switching_plus.out/table_49x100_300K_alpha0.01_100ns_2nm.txt'
    # inFileName = 'biaxial_island_switching_plus.out/table_65x100_350K_alpha0.01_0.5µs_2nm.txt'
    # inFileName = 'biaxial_island_switching_plus.out/table_65x100_350K_alpha0.01_1µs_3.125nm.txt'
    # inFileName = 'biaxial_island_switching_extfield.out/table_50x45_ext0.00015_ext3Pi8_1.25µs_3.125nm.txt'
    # inFileName = 'biaxial_island_switching_extfield.out/table_50x45_ext0_10ns_3.125nm.txt'
    # inFileName = 'biaxial_island_switching_extfield.out/table.txt'
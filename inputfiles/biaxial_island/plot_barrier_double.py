"""
    PLOTS: - 'Energy barrier' AS FUNCTION OF 'Roundness' FOR DIFFERENT 'Ellipse size'
"""
import matplotlib
import os
import numpy as np
import matplotlib.pyplot as plt
import pandas as pd

font = {'size':16}
matplotlib.rc('font', **font)
matplotlib.rcParams['axes.prop_cycle'] = matplotlib.cycler(color=["b", "y", "r", "lightgreen"]) 


def read_mumax3_table(filename):
    """Puts the mumax3 output table in a pandas dataframe"""
    
    table = pd.read_csv(filename, sep='\t')
    table.columns = ' '.join(table.columns).split()[1::2]
    
    return table

def plot(filename, save=False, inset=None, groupBy=None):
    '''
        @param filename [str]: relative path of the mumax table
        @param save [bool] <False>: Whether or not to save the created plot.
        @param inset [list(4 or 6)]: x0, y0, dx, dy, x1, y1, x2, y2 of inset. If not specified, no inset is generated.
            x0, y0, dx, dy specify the physical location and size of the inset on the larger plot, as a ratio.
            x1, y1, x2, y2 specify the range of the plot.
        @param groupBy [str]: If specified, a plot will be made for each separate value in the column named <groupBy>.
    '''
    outDir = 'Figures/Barrier'
    if not os.path.exists(outDir):
        os.makedirs(outDir)
    outFileName = os.path.join(outDir, os.path.splitext(os.path.basename(filename).split('table')[-1])[0]) + '.pdf'

    shape = read_mumax3_table(filename)

    INSET = bool(inset)
    USE_ELECTRONVOLT = True
    USE_ABSOLUTE_VALUE = False

    fig = plt.figure(figsize=(8.0, 5.0))
    ax = fig.add_subplot(111)
    if INSET:
        axins = ax.inset_axes(inset[:4])

    longest_roundnesses = np.array([])
    
    for subset in shape.groupby(groupBy, sort=False):
        size = subset[0]
        subtable = subset[1]
        E_barrier = []
        roundnesses = []
        for subsubset in subtable.groupby("Roundness"):
            roundness = subsubset[0]
            subsubtable = subsubset[1]
            E_demag = subsubtable["E_total"]-subsubtable["E_Zeeman"]
            if USE_ABSOLUTE_VALUE:
                diff = max(E_demag) - min(E_demag)
            elif len(E_demag) == 2:
                diff = E_demag.iloc[1] - E_demag.iloc[0]
            else:
                allowed_energies = E_demag
                allowed_energies = [w for i, w in enumerate(E_demag) if subsubtable["my"].iloc[i] >= -1e-6 and subsubtable["mx"].iloc[i] >= -1e-6]
                # Find index closest to 0°
                index_0 = (np.abs(subsubtable["my"])).idxmin()
                # Find index closest to 45°
                index_45 = (np.abs(np.abs(subsubtable["my"]) - np.sqrt(2)/2)).idxmin()
                # print(E_demag, index_0, index_45)
                # Determine if energy barrier is positive (easy axis 0°) or negative (easy axis 45°)
                sign = 1 if E_demag[index_0] < E_demag[index_45] else -1
                diff = (max(allowed_energies) - min(allowed_energies))*sign
            # print("(%.2f x %.2f nm) Delta E: %.2e J = %.3f eV" % (size, size*roundness, diff, diff/1.602e-19))
            E_barrier.append(diff/1.602e-19 if USE_ELECTRONVOLT else diff)
            roundnesses.append(roundness)
        if len(roundnesses) > len(longest_roundnesses):
            longest_roundnesses = np.array(roundnesses)
        if groupBy == "Size":
            label = '%s nm' % size
        elif groupBy == "Cell_size":
            label = '%s nm' % (size*1e9)
        plt.plot(roundnesses, E_barrier, '.-', label=label)
        if INSET:
            axins.plot(roundnesses, E_barrier, 'o-')

    plt.grid(color='grey', linestyle=':', linewidth=1)
    plt.axhline(0, color='black', linestyle=':', linewidth=1, label=None)
    plt.xlabel(r'Roundness $\rho$')
    plt.ylabel(r'Energy barrier [%s]' % ('eV' if USE_ELECTRONVOLT else 'J'))
    plt.legend()

    if INSET:
        x1, y1, x2, y2 = inset[4:]
        axins.set_xlim(x1, x2)
        axins.set_ylim(y1, y2)
        axins.grid(axis='y', color='grey', linestyle=':', linewidth=1)
        axins.axhline(0, color='black', linestyle=':', linewidth=1, label=None)
        axins.tick_params(labelsize=10)
        axins.set_xticks([i for i in axins.get_xticks() if np.any(np.abs(longest_roundnesses - i) < 1e-6)]) # Only use the major ticks which are a roundness
        axins.set_xticks(longest_roundnesses[np.where(np.logical_and(x1 <= longest_roundnesses, x2 >= longest_roundnesses))[0]], minor=True) # Minor ticks
        ax.indicate_inset_zoom(axins, zorder=10, alpha=0.8)

    plt.gcf().tight_layout()
    if save:
        plt.savefig(outFileName)
    plt.show()

if __name__ == "__main__":
    SAVE = False
    # plot('biaxial_island_shape.out/tablePlus_100_0.1-1_aPi4_B0.001_cell1,2,4nm.txt', save=SAVE, groupBy="Cell_size")
    # plot('biaxial_island_shape.out/tablePlus_32-128_0.1-1_aPi4_B0.001_cell1nm.txt', save=SAVE, groupBy="Size")
    # plot('biaxial_island_shape.out/tablePlus_100_0.1-1_aPi4_B0.001_cell3.125nm.txt', save=SAVE)

    # Those with aPi128 include the second-order anisotropy, at the cost of using a higher field.
    plot('biaxial_island_shape.out/tablePlus_100_0.1-1_aPi128_B0.01_cell1,2,4nm.txt', save=SAVE, groupBy="Cell_size", inset=[0.6, 0.1, 0.36, 0.37, 0.45, -0.15, 0.55, 0.15])
    # plot('biaxial_island_shape.out/tablePlus_32,64,128_0.1-1_aPi128_B0.01_cell1nm.txt', save=SAVE, groupBy="Size")

    # plot('biaxial_island_shape.out/table.txt', save=SAVE)
    
"""
    PLOTS: - 'Energy' AS FUNCTION OF 'Relaxed magnetization angle' FOR DIFFERENT 'External magnetic field strength'
           OR 'Energy' AS FUNCTION OF 'External magnetic field angle' FOR DIFFERENT 'External magnetic field strength'
           - 'Energy barrier' AS FUNCTION OF 'External magnetic field strength'
"""
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

if __name__ == "__main__":
    SAVE1 = False
    SAVE2 = False
    USE_ELECTRONVOLT = True
    USE_RELAXED_ANGLE = True
    
    inFileName = 'biaxial_island_shape_field.out/tablePlus_65_B25-0.001-div4_a128Pi_plotOptimized.txt'
    outDir1 = 'Figures/BarrierLandscape'
    if not os.path.exists(outDir1):
        os.makedirs(outDir1)
    outFileName1 = os.path.join(outDir1, os.path.splitext(os.path.basename(inFileName).split('table')[-1])[0]) + '.pdf'

    shape = read_mumax3_table(inFileName)

    fig = plt.figure(figsize=(8.0, 6.0))
    legend = []
    fields, E_min, E_max = [], [], []
    ELECTRONVOLT_FACTOR = 1.602e-19 if USE_ELECTRONVOLT else 1
    for subset in shape.groupby("Field", sort=False):
        field = subset[0]
        subtable = subset[1]
        E_demag = subtable["E_total"]-subtable["E_Zeeman"]
        E_demag /= ELECTRONVOLT_FACTOR
        diff = max(E_demag) - min(E_demag)
        print("(%.3f T) Delta E: %.2e %s" % (field, diff, 'eV' if USE_ELECTRONVOLT else 'J'))
        legend.append('%.4f T' % field)
        # legend.append(('%.1f T' % (field)) if field > 1 else ('%d mT' % (field*1e3)))

        if USE_RELAXED_ANGLE:
            angles = np.arctan2(subtable["my"], subtable["mx"])*180/np.pi
            if 1: # Plot angles successively (180 -> 181 instead of 180 -> -179)
                previousAngles = np.array(angles)[:-1]
                nextAngles = np.array(angles)[1:]
                offsets = (np.abs(previousAngles - nextAngles) > 180)*((previousAngles > nextAngles)*2 - 1)*360
                offset = np.cumsum(offsets)
                fancyAngles = np.append([angles.iloc[0]], angles.iloc[1:] + offset)
                angles = fancyAngles
            
            angles, E_demag = zip(*sorted(zip(angles, E_demag))) # Sort angles
            plt.scatter(angles, E_demag) # Energy as function of relaxed magnetization angle
        else: 
            plt.scatter(subtable["Angle"], E_demag) # Energy as function of external magnetic field angle
        fields.append(field)
        E_min.append(min(E_demag))
        E_max.append(max(E_demag))
    E_barrier = [E_max[i] - E_min[i] for i, _ in enumerate(fields)]

    plt.grid(color='grey', linestyle=':', linewidth=1)
    plt.legend(legend)
    plt.xlim([0,90])
    plt.ylim([8.1e-19/ELECTRONVOLT_FACTOR, 8.8e-19/ELECTRONVOLT_FACTOR])
    plt.xlabel(r"Relaxed magnetization angle $\theta$ [°]" if USE_RELAXED_ANGLE else r"External magnetic field angle $\theta$ [°]")
    plt.ylabel(r"Energy [%s]" % ('eV' if USE_ELECTRONVOLT else 'J'))
    plt.title(r"65x100 nm double-ellipse")

    # Show plot
    plt.gcf().tight_layout()
    if SAVE1:
        plt.savefig(outFileName1)
    plt.show()


    outDir2 = 'Figures/Barrier'
    if not os.path.exists(outDir2):
        os.makedirs(outDir2)
    outFileName2 = os.path.join(outDir2, 'Field' + os.path.splitext(os.path.basename(inFileName).split('table')[-1])[0]) + '.pdf'

    fig = plt.figure(figsize=(8.0, 5.0))
    # Show trend of E_min and E_max
    plt.plot(fields, E_min)
    plt.scatter(fields, E_min)
    plt.plot(fields, E_max)
    plt.scatter(fields, E_max)
    plt.plot(fields, E_barrier)
    plt.scatter(fields, E_barrier)
    plt.legend([r'$E_{min}$', r'$E_{max}$', r'$E_{barrier}$'])
    plt.xscale("log")
    plt.xlabel(r"$B_{ext}$ [T]")
    plt.ylabel(r"Energy barrier [%s]" % ('eV' if USE_ELECTRONVOLT else 'J'))
    plt.title(r"65x100 nm double-ellipse")
    plt.gcf().tight_layout()
    if SAVE2:
        plt.savefig(outFileName2)
    plt.show()
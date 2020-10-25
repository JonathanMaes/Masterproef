"""
    PLOTS: - 'Energy' AS FUNCTION OF 'Relaxed magnetization angle' FOR DIFFERENT 'External magnetic field strength'
           OR 'Energy' AS FUNCTION OF 'External magnetic field angle' FOR DIFFERENT 'External magnetic field strength'
           - 'Energy barrier' AS FUNCTION OF 'External magnetic field strength'
"""
import numpy as np
import matplotlib.pyplot as plt
import pandas as pd


def read_mumax3_table(filename):
    """Puts the mumax3 output table in a pandas dataframe"""
    
    table = pd.read_csv(filename, sep='\t')
    table.columns = ' '.join(table.columns).split()[1::2]
    
    return table

if __name__ == "__main__":
    # Shape anisotropy: Plus
    shape = read_mumax3_table('biaxial_island_shape_field.out/tablePlus_65_T100-0.001_a512Pi.txt')
    legend = []
    fields, E_min, E_max = [], [], []
    for subset in shape.groupby("Field"):
        field = subset[0]
        subtable = subset[1]
        E_demag = subtable["E_total"]-subtable["E_Zeeman"]
        diff = max(E_demag) - min(E_demag)
        print("(%.3f T) Delta E: %.2e J = %.3f eV" % (field, diff, diff/1.602e-19))
        legend.append('%.4f T' % field)

        angles = np.arctan2(subtable["my"], subtable["mx"])*180/np.pi
        if 1: # Plot angles successively (180 -> 181 instead of 180 -> -179)
            previousAngles = np.array(angles)[:-1]
            nextAngles = np.array(angles)[1:]
            offsets = (np.abs(previousAngles - nextAngles) > 180)*((previousAngles > nextAngles)*2 - 1)*360
            offset = np.cumsum(offsets)
            fancyAngles = np.append([angles.iloc[0]], angles.iloc[1:] + offset)
            angles = fancyAngles
        
        angles, E_demag = zip(*sorted(zip(angles, E_demag))) # Sort angles
        plt.plot(angles, E_demag) # Energy as function of relaxed magnetization angle
        # plt.plot(subtable["Angle"], E_demag) # Energy as function of external magnetic field angle
        fields.append(field)
        E_min.append(min(E_demag))
        E_max.append(max(E_demag))
    E_barrier = [E_max[i] - E_min[i] for i, _ in enumerate(fields)]
    plt.legend(legend)
    plt.xlabel(r"Relaxed magnetization angle $\theta$ [°]")
    # plt.xlabel(r"External magnetic field angle $\theta$ [°]")
    plt.ylabel(r"Energy [J]")
    plt.title(r"65x100 nm double-ellipse")

    # Show plot
    plt.show()

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
    plt.ylabel(r"Energy barrier [J]")
    plt.title(r"65x100 nm double-ellipse")
    plt.show()
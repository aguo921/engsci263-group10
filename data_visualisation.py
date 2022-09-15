import matplotlib.pyplot as plt
from mpl_toolkits import axisartist
from mpl_toolkits.axes_grid1 import host_subplot
import numpy as np

def plot_data():
    """ Plot the given pressure, subsidence and mass extraction data. """

    # extract data from given files
    tP, P = np.genfromtxt("sb_pres.txt", delimiter=",", skip_header=1).T
    tU, U = np.genfromtxt("sb_disp.txt", delimiter=",", skip_header=1).T
    tq, q = np.genfromtxt("sb_mass.txt", delimiter=",", skip_header=1).T

    # create axes
    host = host_subplot(111,axes_class=axisartist.Axes)
    plt.subplots_adjust(right=0.75)

    par1 = host.twinx()
    par2 = host.twinx()

    par2.axis["right"] = par2.new_fixed_axis(loc="right",offset=(60,0))
    par1.axis["right"].toggle(all=True)
    par2.axis["right"].toggle(all=True)

    # plot the data
    host.plot(tU,U,label="subsidence")
    par1.plot(tP,P,label="pressure")
    par2.plot(tq,q,label="mass extraction")

    # set axes range
    host.set_xlim(1952, 2014)
    host.set_ylim(0, 16)
    par1.set_ylim(20, 60)
    par2.set_ylim(0, 1800)

    # set axes labels
    host.set_xlabel("time [years]")
    host.set_ylabel("subsidence [m]")
    par1.set_ylabel("pressure [bars]")
    par2.set_ylabel("mass extraction [kg/s]")

    # set title
    host.set_title("Historical data from Wairakei Geothermal system")

    # show legend
    host.legend(loc="lower right")

    # display plots
    plt.show()


if __name__ == "__main__":
    plot_data()
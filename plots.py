import numpy as np
import matplotlib.pyplot as plt
import matplotlib.ticker as mticker

def plot():
    x, y, z = np.loadtxt('testy/wykres.txt', unpack=True, delimiter=',')
    fig = plt.figure()

    ax = plt.axes(projection='3d')
    my_cmap = plt.get_cmap('hot')

    trisurf = ax.plot_trisurf(x, y, z, cmap = my_cmap,
                         linewidth = 0.2,
                         antialiased = True,
                         edgecolor = 'grey')
    #fig.colorbar(trisurf, ax = ax, shrink = 0.5, aspect = 5)
   
    ax.set_title('gr120.atsp\n')
    ax.set_xlabel('probability of mutation [%]')
    ax.set_ylabel('iterations without update')
    ax.set_zlabel('time [s]')
    #ax.set_xlabel('size of population')
    #ax.set_ylabel('iterations without update')
    #ax.set_ylabel('probability of muatation [%]')
    #ax.set_ylabel('iterations without update')

    formatter = mticker.ScalarFormatter()
    ax.xaxis.set_major_formatter(formatter)
    ax.xaxis.set_major_locator(mticker.FixedLocator(x))
    ax.yaxis.set_major_locator(mticker.FixedLocator(y))

    plt.show()


if __name__ == '__main__':
    plot()
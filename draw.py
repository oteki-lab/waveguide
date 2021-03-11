from numpy import pi, sqrt, exp, sin, cos, arctan
import numpy as np
import matplotlib.pyplot as plt

def graph1(L0, wl0, b, Gamma):
    fig = plt.figure()
    ax1 = fig.add_subplot(1, 2, 1)
    for k, _ in enumerate(wl0):
        ax1.plot(L0*1e9, b.transpose()[k], label='{} um'.format(int(wl0[k]*1e6)))
    plt.xlim(0, max(L0)*1e9)
    plt.ylim(0,1)
    plt.legend(bbox_to_anchor=(0, 1), loc='upper left', borderaxespad=0, fontsize=9, prop={'size':8,}, frameon=False, title='wavelength')
    ax2 = fig.add_subplot(1, 2, 2)
    for k, _ in enumerate(wl0):
        ax2.plot(L0*1e9, Gamma.transpose()[k], label='{} um'.format(int(wl0[k]*1e6)))
    plt.xlim(0, max(L0)*1e9)
    plt.ylim(0,1)
    plt.legend(bbox_to_anchor=(0, 1), loc='upper left', borderaxespad=0, fontsize=9, prop={'size':8,}, frameon=False, title='wavelength')
    plt.show()

def graph2(L0, wl0, j, E, y):
    ylim = [0, 1.1]
    fig = plt.figure()
    ax3 = fig.add_subplot(1, 1, 1)
    for k, _ in enumerate(wl0):
        ax3.plot(y[j]*1e9, E[j][k], label='{} um'.format(int(wl0[k]*1e6)))
    ax3.plot(np.array([-L0[j]/2, -L0[j]/2])*1e9, ylim, color='black', linestyle='dashed')
    ax3.plot(np.array([L0[j]/2, L0[j]/2])*1e9, ylim, color='black', linestyle='dashed')
    plt.xlim(min(y[j])*1e9, max(y[j])*1e9)
    plt.ylim(ylim)
    #plt.yscale('log')
    plt.legend(borderaxespad=0, fontsize=9, prop={'size':8,}, frameon=False, title='wavelength')
    plt.show()

from numpy import pi, sqrt, exp, sin, cos, arctan
import numpy as np
import matplotlib.pyplot as plt

def graph1(L0, wl0, b, Gamma):
    fig = plt.figure()
    ax1 = fig.add_subplot(121, xlabel='Core width (nm)', ylabel='Normalized guide index')
    for k, _ in enumerate(wl0):
#        ax1.plot(L0*1e9, b.transpose()[k], label='{:.2f} um'.format(wl0[k]*1e6))
        ax1.plot(L0*1e9, b.transpose()[k], label='{}'.format(wl0[k]))
    plt.xlim(0, max(L0)*1e9)
    plt.ylim(0,1)
    ax1.set_yticks([0.0, 0.1, 0.2, 0.3, 0.4, 0.5 ,0.6 ,0.7 ,0.8, 0.9, 1.0])
    ax1.grid(which = "major", axis = "y", color = "gray", alpha = 0.8, linestyle = "--", linewidth = 1)
#    plt.legend(bbox_to_anchor=(1, 0), loc='lower right', borderaxespad=0, fontsize=9, prop={'size':8,}, frameon=False, title='wavelength')
    plt.legend(bbox_to_anchor=(1, 0), loc='lower right', borderaxespad=0, fontsize=9, prop={'size':8,}, frameon=False, title='mode')
    ax2 = fig.add_subplot(122, xlabel='Core width (nm)', ylabel='Confinement factor')
    for k, _ in enumerate(wl0):
#        ax2.plot(L0*1e9, Gamma.transpose()[k], label='{:.2f} um'.format(wl0[k]*1e6))
        ax2.plot(L0*1e9, Gamma.transpose()[k], label='{}'.format(wl0[k]))
    plt.xlim(0, max(L0)*1e9)
    plt.ylim(0,1)
    ax2.set_yticks([0.0, 0.1, 0.2, 0.3, 0.4, 0.5 ,0.6 ,0.7 ,0.8, 0.9, 1.0])
    ax2.grid(which = "major", axis = "y", color = "gray", alpha = 0.8, linestyle = "--", linewidth = 1)
#    plt.legend(bbox_to_anchor=(1, 0), loc='lower right', borderaxespad=0, fontsize=9, prop={'size':8,}, frameon=False, title='wavelength')
    plt.legend(bbox_to_anchor=(1, 0), loc='lower right', borderaxespad=0, fontsize=9, prop={'size':8,}, frameon=False, title='mode')
    fig.tight_layout()
    plt.show()

def graph2(L0, wl0, j, E, y):
    ylim = [0, 1.1]
    fig = plt.figure()
    ax3 = fig.add_subplot(111, xlabel='y (nm)', ylabel='Ey')
    for k, _ in enumerate(wl0):
#        ax3.plot(y[j]*1e9, E[j][k], label='{:.2f} um'.format(wl0[k]*1e6))
        ax3.plot(y[j]*1e9, E[j][k], label='{}'.format(wl0[k]))
    ax3.plot(np.array([-L0[j]/2, -L0[j]/2])*1e9, ylim, color='black', linestyle='dashed')
    ax3.plot(np.array([L0[j]/2, L0[j]/2])*1e9, ylim, color='black', linestyle='dashed')
    plt.xlim(min(y[j])*1e9, max(y[j])*1e9)
    plt.ylim(ylim)
    #plt.yscale('log')
#    plt.legend(borderaxespad=0, fontsize=9, prop={'size':8,}, frameon=False, title='wavelength')
    plt.legend(borderaxespad=0, fontsize=9, prop={'size':8,}, frameon=False, title='mode')
    plt.show()

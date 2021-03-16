from numpy import pi, sqrt, exp, sin, cos, arctan
import numpy as np
import matplotlib.pyplot as plt


def graph1(m, w, d, b, Gamma, mode):
    fig = plt.figure()
    ax1 = fig.add_subplot(121, xlabel='Core width (nm)', ylabel='Normalized guide index')
    if type(m) is np.ndarray:
        for i, _ in enumerate(m):
            ax1.plot(d*1e9, b[i][w], label='{}'.format(m[i]))
            title='mode'
    else:
        for j, _ in enumerate(w):
            ax1.plot(d*1e9, b[m][j], label='{:.2f} um'.format(w[j]*1e6))
            title='wavelength'
    plt.xlim(0, max(d)*1e9)
    plt.ylim(0,1)
    ax1.set_yticks([0.0, 0.1, 0.2, 0.3, 0.4, 0.5 ,0.6 ,0.7 ,0.8, 0.9, 1.0])
    ax1.grid(which = "major", axis = "y", color = "gray", alpha = 0.8, linestyle = "--", linewidth = 1)
    plt.legend(bbox_to_anchor=(1, 0), loc='lower right', borderaxespad=0, fontsize=9, prop={'size':8,}, frameon=False, title=title)
    #ax1.text(0, 0.95, mode)
    ax2 = fig.add_subplot(122, xlabel='Core width (nm)', ylabel='Confinement factor')
    if type(m) is np.ndarray:
        for i, _ in enumerate(m):
            ax2.plot(d*1e9, Gamma[i][w], label='{}'.format(m[i]))
            title='mode'
    else:
        for j, _ in enumerate(w):
            ax2.plot(d*1e9, Gamma[m][j], label='{:.2f} um'.format(w[j]*1e6))
            title='wavelength'
    plt.xlim(0, max(d)*1e9)
    plt.ylim(0,1)
    ax2.set_yticks([0.0, 0.1, 0.2, 0.3, 0.4, 0.5 ,0.6 ,0.7 ,0.8, 0.9, 1.0])
    ax2.grid(which = "major", axis = "y", color = "gray", alpha = 0.8, linestyle = "--", linewidth = 1)
    plt.legend(bbox_to_anchor=(1, 0), loc='lower right', borderaxespad=0, fontsize=9, prop={'size':8,}, frameon=False, title=title)
    #ax2.text(0, 0.95, mode)
    fig.tight_layout()
    plt.show()

def graph2(m, w, y, E, boundaries, mode):
    fig = plt.figure()
    ax = fig.add_subplot(111, xlabel='y (nm)', ylabel='Ey')
    if type(m) is np.ndarray:
        for i, _ in enumerate(m):
            ax.plot(y*1e9, E[i][w], label='{}'.format(m[i]))
            title='mode'
    else:
        for j, _ in enumerate(w):
            ax.plot(y*1e9, E[m][j], label='{:.2f} um'.format(w[j]*1e6))
            title='wavelength'
    ylim = [-1.1, 1.1]
    for boundary in boundaries:
        ax.plot([boundary, boundary], ylim, color='black', linestyle='dashed')
    plt.xlim(np.array([min(y), max(y)])*1e9)
    plt.ylim(ylim)
    plt.legend(borderaxespad=0, fontsize=9, prop={'size':8,}, frameon=False, title=title)
    #ax.text(0, 1.0, mode)
    plt.show()

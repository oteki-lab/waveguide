from numpy import pi, sqrt, exp, sin, cos, arctan
import numpy as np
import matplotlib.pyplot as plt
import n

wl = 1100e-9
L = np.arange(10, 1800+10, 10)*1e-9
Nt = n.AlGaAs(0.8, wl)
Ng = n.AlGaAs(0.0, wl)
Ns = n.AlGaAs(0.8, wl)
m=0


def left(x):
    return 2*pi*L/wl*sqrt(Ng**2-x**2)
def right(x):
    return arctan(sqrt((x**2-Nt**2)/(Ng**2-x**2))) + arctan(sqrt((x**2-Ns**2)/(Ng**2-x**2))) + m*pi
def h(x):
    return left(x)-right(x)


def BinarySearch(eq, upper,lower):
    err = 1e-10
    up=upper-err*np.ones(len(L))
    low=lower+err*np.ones(len(L))
    limit=10000
    for count in range(limit):
        mid=(up+low)/2  #数値解
        y=eq(mid)       #関数値
        if all(abs(y)<err) or count>limit:      #解が発見された
            break
        flag = eq(low)*y<0
        up[flag]=mid[flag]    #解は下限と中間点の間にある
        low[~flag]=mid[~flag] #解は上限と中間点の間にある
    return mid

k0 = 2*pi/wl
Ne = BinarySearch(h, Ng, Ns)
E1 = 1/(2*k0*sqrt(Ne**2-Ns**2))
E2 = L/2 + sin(2*k0*sqrt(Ng**2-Ne**2)*L)/(4*k0*sqrt(Ng**2-Ne**2))
E3 = (cos(k0*sqrt(Ng**2-Ne**2)*L))**2/(2*k0*sqrt(Ne**2-Nt**2))
Gamma = E2/(E1+E2+E3)
Gamma[Ng-Ne<1e-5] = 0.0
b = (Ne**2-Nt**2)/(Ng**2-Nt**2)
b[Ng-Ne<1e-5] = 0.0
As = k0*sqrt(Ne**2-Ns**2)
At = k0*sqrt(Ne**2-Nt**2)
phi = arctan(-As/(k0*sqrt(Ng**2-Ne**2)))

fig = plt.figure()
ax1 = fig.add_subplot(1, 2, 1)
ax1.plot(L*1e9, b)
plt.xlim(0, max(L)*1e9)
plt.ylim(0,1)
ax2 = fig.add_subplot(1, 2, 2)
ax2.plot(L*1e9, Gamma)
plt.xlim(0, max(L)*1e9)
plt.ylim(0,1)
plt.show()

def f(yi, j):
    if yi<0:
        return cos(phi[j])*exp(As[j]*yi)
    elif yi<L[j]:
        return cos(k0*sqrt(Ng**2-Ne[j]**2)*yi+phi[j])
    else:
        return cos(k0*sqrt(Ng**2-Ne[j]**2)*L[j]+phi[j])*exp(-At[j]*(yi-L[j]))

j = -1
y = np.arange(-1000, L[j]*1e9+1010, 10)*1e-9
x = y-L[j]/2
E = np.array([f(yi, j) for yi in y])

ylim = [0, 1.1]
fig = plt.figure()
ax3 = fig.add_subplot(1, 1, 1)
ax3.plot(x*1e9, E)
ax3.plot(np.array([-L[j]/2, -L[j]/2])*1e9, ylim, color='black', linestyle='dashed')
ax3.plot(np.array([L[j]/2, L[j]/2])*1e9, ylim, color='black', linestyle='dashed')
plt.xlim(min(x)*1e9, max(x)*1e9)
plt.ylim(ylim)
#plt.yscale('log')
plt.show()

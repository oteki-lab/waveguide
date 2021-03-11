from numpy import pi, sqrt, exp, sin, cos, arctan
import numpy as np
import matplotlib.pyplot as plt
import n
import draw

wl0 = np.arange(1000, 5000+1000, 1000)*1e-9    #1200e-9
L0 = np.arange(10, 2000+10, 10)*1e-9
wl = np.tile(wl0, (len(L0), 1))
L = np.tile(L0, (len(wl0), 1)).transpose()
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
    up=upper-err
    low=lower+err
    limit=10000
    for count in range(limit):
        mid=(up+low)/2  #数値解
        y=eq(mid)       #関数値
        if np.all(abs(y)<err*np.ones(y.shape)) or count>limit:      #解が発見された
            break
        flag = eq(low)*y<0
        up[flag]=mid[flag]    #解は下限と中間点の間にある
        low[~flag]=mid[~flag] #解は上限と中間点の間にある
    return mid


Ne = BinarySearch(h, Ng, Ns)
k0 = 2*pi/wl
ky = k0*sqrt(Ng**2-Ne**2)
As = k0*sqrt(Ne**2-Ns**2)
At = k0*sqrt(Ne**2-Nt**2)
phi = arctan(-As/(ky))
E1 = 1/(2*As)
E2 = L/2 + sin(2*ky*L)/(4*ky)
E3 = (cos(ky*L))**2/(2*At)
Gamma = E2/(E1+E2+E3)
Gamma[Ng-Ne<1e-5] = 0.0
b = (Ne**2-Nt**2)/(Ng**2-Nt**2)
b[Ng-Ne<1e-5] = 0.0

draw.graph1(L0, wl0, b, Gamma)

def f(yi, j, k):
    if yi<0:
        return cos(phi[j][k])*exp(As[j][k]*yi)
    elif yi<L[j][k]:
        return cos(k0[j][k]*sqrt(Ng[j][k]**2-Ne[j][k]**2)*yi+phi[j][k])
    else:
        return cos(k0[j][k]*sqrt(Ng[j][k]**2-Ne[j][k]**2)*L[j][k]+phi[j][k])*exp(-At[j][k]*(yi-L[j][k]))

div = 1000
y = np.empty((len(L0), div))
E = np.empty((len(L0), len(wl0), div))
for j, Lj in enumerate(L0):
    for k, wlk in enumerate(wl0):
        y[j] = np.linspace(-L0[j]*3/2, L0[j]*3/2, div)
        E[j][k] = np.array([f(yi+L0[j]/2, j, k) for yi in y[j]])

j = -1
draw.graph2(L0, wl0, j, E, y)
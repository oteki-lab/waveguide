from numpy import pi, sqrt, arctan
import numpy as np

wl = 1100e-9
L = np.array([100e-9, 1000e-9])
Nt = 3.3
Ng = 3.5
Ns = 3.3
m=0


def f(x):
    return 2*pi*L/wl*sqrt(Ng**2-x**2)
def g(x):
    return arctan(sqrt((x**2-Nt**2)/(Ng**2-x**2))) + arctan(sqrt((x**2-Ns**2)/(Ng**2-x**2))) + m*pi
def h(x):
    return f(x)-g(x)


def BinarySearch(eq, upper,lower):
    err = 1e-8
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

Neff = BinarySearch(h, Ng, Ns)

print(Neff)

from numpy import pi, sqrt, exp, sin, cos, arctan
import numpy as np
import matplotlib.pyplot as plt
import n
import draw

mode = 'TE'
m = np.arange(0, 5+1)
w = np.arange(1000, 5000+500, 500)*1e-9
d = np.arange(10, 2000+10, 10)*1e-9

M = np.array([np.tile(m, (len(w),1))]*len(d)).transpose()
D = np.array([np.tile(d, (len(w),1))]*len(m))
W = np.array([np.tile(w, (len(m),1)).transpose()]*len(d)).transpose()

Nt = n.AlGaAs(0.8, W)
Ng = n.AlGaAs(0.0, W)
Ns = n.AlGaAs(0.8, W)

def left(x):
    return 2*pi*D/W*sqrt(Ng**2-x**2)
def right(x, mode):
    if mode == 'TE':
        return arctan(sqrt((x**2-Nt**2)/(Ng**2-x**2))) + arctan(sqrt((x**2-Ns**2)/(Ng**2-x**2))) + M*pi
    else:
        return arctan((Ng/Nt)**2*sqrt((x**2-Nt**2)/(Ng**2-x**2))) + arctan(sqrt((Ng/Ns)**2*(x**2-Ns**2)/(Ng**2-x**2))) + M*pi
def h(x, mode):
    return left(x)-right(x, mode)

def BinarySearch(eq, upper, lower, mode):
    err = 1e-10
    up=upper-err
    low=lower+err
    limit=10000
    for count in range(limit):
        mid=(up+low)/2  #数値解
        y=eq(mid, mode)       #関数値
        if np.all(abs(y)<err*np.ones(y.shape)) or count>limit:      #解が発見された
            break
        flag = eq(low, mode)*y<0
        up[flag]=mid[flag]    #解は下限と中間点の間にある
        low[~flag]=mid[~flag] #解は上限と中間点の間にある
    return mid

Ne = BinarySearch(h, Ng, Ns, mode)

k0 = 2*pi/W
ky = k0*sqrt(Ng**2-Ne**2)
As = k0*sqrt(Ne**2-Ns**2)
At = k0*sqrt(Ne**2-Nt**2)
phi = arctan(-As/ky)

E1, E2, E3 = [1/(2*As), D/2+sin(2*ky*D)/(4*ky), (cos(ky*D))**2/(2*At)]
Gamma = E2/(E1+E2+E3)
Gamma[Ng-Ne<1e-5] = 0.0
b = (Ne**2-Nt**2)/(Ng**2-Nt**2)
b[Ng-Ne<1e-5] = 0.0

i = 0
j = 0

draw.graph1(m, j, d, b, Gamma, mode)


def f(yl, i, j, k):
    if yl<0:
        return cos(phi[i][j][k])*exp(As[i][j][k]*yl)
    elif yl<D[i][j][k]:
        return cos(k0[i][j][k]*sqrt(Ng[i][j][k]**2-Ne[i][j][k]**2)*yl+phi[i][j][k])
    else:
        return cos(k0[i][j][k]*sqrt(Ng[i][j][k]**2-Ne[i][j][k]**2)*D[i][j][k]+phi[i][j][k])*exp(-At[i][j][k]*(yl-D[i][j][k]))

k = -1
start = -d[k]*3/2
end = d[k]*3/2
div = 1000
y = np.linspace(start, end, div)

E = np.array([[[
    f(yl+d[k]/2, i, j, k) 
    for yl in y] 
    for j, _ in enumerate(w)] 
    for i, _ in enumerate(m)]
)

i = -1
j = 0
boundaries = np.array([-d[k]/2, d[k]/2])*1e9

draw.graph2(m, j, y, E, boundaries, mode)

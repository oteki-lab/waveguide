from numpy import pi, sqrt, exp, sin, cos, arctan
import numpy as np

def AlGaAs(x, l):
    h = 6.626e-34   # [Js]
    c = 2.998e8     # [m/s]
    e = 1.60218e-19 # [J] = 1 eV

    A0 = 6.3+19.0*x
    B0 = 9.4-10.2*x
    E0 = 1.425 + 1.155*x + 0.37*x**2    # [eV]
    Ed0 = 1.765 + 1.115*x + 0.37*x**2    # [eV]

    chi = h*c/(l*E0)/e
    chi_s0 = h*c/(l*Ed0)/e

    def f(y):
        return (2-sqrt(1+y)-sqrt(1-y))/y**2

    n = sqrt( A0*(f(chi)+f(chi_s0)/2*(E0/Ed0)**(3/2)) + B0 )

    return n
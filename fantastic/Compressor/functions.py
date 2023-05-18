import numpy as np
from numpy.polynomial import Polynomial
import globals as gb

def dh0(delta_p, rho, eff):
    enthalpy = delta_p/(rho * eff)
    return enthalpy

def v_radial(mdot, r, rho, h):
    v_r = mdot/(2 * np.pi * r * h * rho)
    return v_r

def euler(dh0, omeg, r1, r2, vtht1=0):
    vtht2 = (dh0/(omeg) - (r1*vtht1))/r2
    return vtht2

def omega_opt(beta_2, V_r1, r_2, sigma, dh0):
    omeg = Polynomial([-dh0, V_r1 * (r_2 ** 2 )* np.tan(beta_2), sigma * r_2 ** 2]).roots()[1]
    return omeg

import numpy as np
from numpy.polynomial import Polynomial
import globals as gb
from matplotlib import pyplot as plt

def dh0(delta_p, rho, eff):
    enthalpy = delta_p/(rho * eff)
    return enthalpy

def v_radial(mdot, r, rho, h):
    v_r = mdot/(2 * np.pi * r * h * rho)
    return v_r

def euler(dh0, omeg, r1, r2, vtht1=0):
    vtht2 = (dh0/(omeg) - (r1*vtht1))/r2
    return vtht2

def omega_opt(beta_2, V_r1, r_2, sigma, dh0, omega):
    opt = Polynomial([-dh0, V_r1 * (r_2 ** 2 )* np.tan(beta_2), sigma * r_2 ** 2]).roots()[1] - omega
    return opt

def beta_mid(beta_1, beta_2):
    beta_mid = np.arctan(0.5*(np.arctan(beta_1) + np.arctan(beta_2)))
    return beta_mid

def rotor_plot(r_1, r_2, beta1, beta2, step):
    r_val = np.arange(r_1, r_2, step)
    theta = np.zeros(len(r_val))
    c1 = (np.tan(beta1) - np.tan(beta2))/(r_1**2/2 - r_2**2/2)
    c2 = np.tan(beta2) - c1 * r_2 ** 2/2
    b = np.arctan(c1/2 * r_val ** 2 + c2)
    for j in range(1, len(r_val)):
        theta[j] = theta[j-1] + np.tan(b[j])/r_val[j] * step
    plt.polar(theta, r_val, 'black')
    plt.yticks([])
    plt.ylim(0, None)
    plt.show()
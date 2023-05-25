import numpy as np
from numpy.polynomial import Polynomial
import globals as gb
from matplotlib import pyplot as plt
import cmath
import numpy as np
import pandas as pd
from numpy.polynomial import Polynomial
from scipy.integrate import quad
from scipy.interpolate import UnivariateSpline

def dh0(delta_p, rho, eff):
    enthalpy = delta_p/(rho * eff)
    return enthalpy

def v_radial(mdot, r, rho, h):
    v_r = mdot/(2 * np.pi * r * h * rho)
    return v_r

def omega_opt(T_loss_0, eta_c, eta_t, vacuum_power_0, suction_box_rise, dp_guess_comp,  V2_r, beta_2, r_2, sigma0, rho):
    count2 = 0
    dh0_c = dh0(dp_guess_comp,rho,eta_c)
    power = vacuum_power_0
    T_loss = T_loss_0
    while count2 < 30:
        a = sigma0 * (r_2**2)
        b = V2_r * (r_2 ** 2) * np.tan(beta_2)
        c = -dh0_c
        d = (b**2) - (4*a*c)
        sol2 = (-b+cmath.sqrt(d))/(2*a)
        omega = np.real(sol2)
        w_turbine, eta_m, eta_ov = omega_calc(T_loss, eta_c, eta_t, power, omega)
        dp_new_comp = suction_box_rise*eta_ov/(1-eta_ov)
        dh0_c = dh0(dp_new_comp,rho,eta_c)
        T_loss = T_data(omega)
        count2 += 1
    return omega

def T_analytic(omega, rho = gb.rho, r_2 = gb.r_2, f_0 = 0.85, nu=100e-6, d_m=0.025):
        C_f = 0.0025
        n = omega * 60/(2 * np.pi)
        M_0 = f_0 * 1e3 * (nu * n) ** (2 / 3) * d_m ** 3
        T_1 = 0.2 * rho *np.pi * C_f * ((r_2)**5) * (omega**2) 
        T_ov = 2 * M_0 + 4 * T_1
        return T_ov

def T_data(omega, I=0.002456):
    t = (1.20979252149514 * (omega**3)/(10**7))-(0.000331891054100759*(omega**2))+(0.332356972914549*omega)-72.7879697323764
    domega_dt = (0.0260*t**2)-(0.8032*t)+14.2316
    T_loss = I*domega_dt
    return T_loss    

def omega_calc(T_loss, eta_c, eta_t, vacuum_power, omega_c):
    P_loss = T_loss * omega_c
    eta_m = (vacuum_power*eta_t - P_loss)/(eta_t*(vacuum_power - (eta_c*P_loss)))
    eta_ov = eta_t * eta_c * eta_m
    w_turbine = vacuum_power * eta_t/(1-eta_ov)
    return w_turbine, eta_m, eta_ov

def velocity2_triangles(r1, r2, v1_r, v2_r, beta2, omega, Nb=0):
    deg = np.pi/180
    if Nb == 0:
        sigma = 0.85
    else:
        sigma = 1 - (np.sqrt(np.cos(beta2))/(Nb**(0.7)))
    v2_thet = (omega*r2*sigma) + v2_r*np.tan(beta2)

    Ub_1 = omega * r1
    Ub_2 = omega * r2
    # Calculate relative velocities
    v1_rel = np.sqrt(v1_r ** 2 + Ub_1 ** 2)
    v2_rel = np.sqrt(v2_r ** 2 + (v2_thet - Ub_2) ** 2)
    if np.arctan(-Ub_1 / v1_r) >= 0:
        beta1 = np.arctan(-Ub_1 / v1_r) - gb.inlet_angle_rotor * deg
    else:
        beta1 = np.arctan(-Ub_1 / v1_r) + gb.inlet_angle_rotor * deg

    #De Haller test
    if v2_rel/v1_rel < 1/3:
        text = 'Invalid selection'
        print(text)
    return(beta1, v2_thet, v1_rel, v2_rel)

def blades_plot(r_1, r_2, beta1, beta2, step):
    r_val = np.arange(r_1, r_2, step)
    theta = np.zeros(len(r_val))
    c1 = (np.tan(beta1) - np.tan(beta2))/(r_1**2/2 - r_2**2/2)
    c2 = np.tan(beta2) - c1/2 * r_2 ** 2
    b = np.arctan(c1/2 * r_val ** 2 + c2)
    for j in range(1, len(r_val)):
        theta[j] = theta[j-1] + np.tan(b[j])/r_val[j] * step
    return(r_val, theta)   

def throat_dist(N_diff, thetadiff, rdiff, blade_h, r_ext, min=0):
    #Computes minimum for min = 1, max otherwise
    thetadiff_2 = thetadiff + (2*np.pi/N_diff)
    r_1 = rdiff * [np.cos(thetadiff), np.sin(thetadiff)]
    r_2 = rdiff * [np.cos(thetadiff_2), np.sin(thetadiff_2)]
    if min == 0:
        Distance = 1E5
        for i in range(len(rdiff)):
            for j in range(len(rdiff)):
                dist2 = ((r_2[0][i] - r_1[0][j])**2) + ((r_2[1][i] - r_1[1][j])**2)
                dist = np.sqrt(dist2)
                if dist < Distance:
                    Distance = dist
                    index = [i, j]
    else:
        Distance = 1E5
        for j in range(len(rdiff)):
            dist2 = ((r_2[0][j] - r_1[0][-1])**2) + ((r_2[1][j] - r_1[1][-1])**2)
            dist = np.sqrt(dist2)
            if dist < Distance:
                Distance = dist
                index = [j, -1]
    area = blade_h * Distance
    return(area, index)

def blade_length(r_array, theta_array, index_start):
    xs = r_array * np.cos(theta_array)
    ys = r_array * np.sin(theta_array)
    blade_length = 0
    for i in np.arange(index_start,len(r_array)-1):
        dr = np.sqrt((xs[i+1] - xs[i])**2 + (ys[i+1] - ys[i])**2)
        blade_length += dr
    return blade_length

def stator_blades(r_3, r_4, beta_3, beta_4, step):
    # plots diffuser blade shapes
    r_val = np.arange(r_3, r_4, step)
    theta = np.zeros(len(r_val))
    b = beta_3 + (beta_4 - beta_3) * (r_val - r_3)/(r_4 - r_3)
    for j in range(1, len(r_val)):
        theta[j] = theta[j-1] + np.tan(b[j])/r_val[j] * step
    return(r_val, theta) 

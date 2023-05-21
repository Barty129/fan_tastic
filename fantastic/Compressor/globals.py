import numpy as np  

m_dot = 0.061 #kg/s
rho = 1.2046 #kg/m3

r_1 = 0.04
r_2 = 0.10 #m
r_4 = 0.15

suction = 11.85 * 10**3 #Pa
suction_power_start = 598.99 #W
delta_p_comp_start = 6.41 * 10**3 #Pa
blade_height = 12 * 10**(-3)#m

v1_theta = 0
eff_comp = 0.6 #initial compressr efficiency estimate
eff_turb = 0.65
sigma_0 = 0.85
G_val = 1.05
T_loss = 0.07168

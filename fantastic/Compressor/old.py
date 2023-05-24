import numpy as np
from numpy.polynomial import Polynomial
import globals as gb
from matplotlib import pyplot as plt
import cmath

def stator_blades(r_3, r_4, beta_3, beta_4, step):
    # plots diffuser blade shapes
    r_val = np.arange(r_3, r_4, step)
    theta = np.zeros(len(r_val))
    b = beta_3 + (beta_4 - beta_3) * (r_val - r_3)/(r_4 - r_3)
    for j in range(1, len(r_val)):
        theta[j] = theta[j-1] + np.tan(b[j])/r_val[j] * step
    return(r_val, theta) 

''''
    file = open("Rot9.txt","w+")
    theta_plot = theta_rot + (8*2*np.pi)/int(Nb)
    for j in range(len(r_rot)):
        x = r_rot[j] * np.cos(theta_plot[j])
        y = r_rot[j] * np.sin(theta_plot[j])
        z = 0
        file.write(str(x))
        file.write(',')
        file.write(str(y))
        file.write(',')
        file.write(str(z))
        file.write('\n')
    file.close()
'''
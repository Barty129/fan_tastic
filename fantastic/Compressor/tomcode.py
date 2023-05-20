from re import L
import numpy as np
import matplotlib.pyplot as plt

# design parameters to be altered
beta2 = 1.2 # trailing edge blade angle in radians
r1 = 0.025 # leading edge radius
inletangle = 5 # angle of incidence of leading edge in degrees
vaneless = 1.1 # ratio of diffuser leading edge to rotor trailing edge radii
arearatio = 3
aspectratio = 10
# numerical parameters
step = 0.001
r4 = 0.15

# design point values
nc = 0.6
h = 0.012
dp0 = 6410
density = 1.2046
mdot = 0.061
r2 = 0.1
dh0 = 1/nc * dp0/density

def bladeplot(N, r1, r2, step, beta1, beta2):
    # plots blade shapes
    r = np.arange(r1, r2, 0.001)
    theta = [[0 for j in range(len(r))] for i in range(N)]
    fig, ax = plt.subplots(subplot_kw={'projection': 'polar'})
    beta = [0 for i in range(len(r))]
    c = (np.tan(beta1) - np.tan(beta2))/(r1**2/2 - r2**2/2)
    d = np.tan(beta2) - c * r2 ** 2/2
    for j in range(len(r)):
        beta[j] = 1*np.arctan(c*r[j]**2/2 + d)
    for i in range(0, N): 
        theta[i][0] = i * 2 * np.pi/N
        for j in range(1, len(r)):
            theta[i][j] = theta[i][j-1] + np.tan(beta[j])/r[j] * step
        ax.plot(theta[i], r)
    ax.grid(True)
    plt.show()
    return(r, theta)

def velocitycalc(N = 0):
    # calculates velocity triangles
    if N != 0:
        sigma = 1 - (np.cos(beta2)**0.5)/N**0.7
    else:
        sigma = 0.85
    Vr1 = mdot/(2 * np.pi * density * r1 * h)
    Vr2 = mdot/(2 * np.pi * density * r2 * h)

    a = sigma * r2 ** 2
    b = r2 * Vr2 * np.tan(beta2)
    c = -dh0
    omega = (-b + (b**2 - 4*a*c)**0.5)/(2*a)

    U2 = omega * r2
    U1 = omega * r1
    Vtheta1 = 0
    VthetaRel1 = -1 * U1
    VthetaRel2ideal = Vr2 * np.tan(beta2)
    Vtheta2ideal = VthetaRel2ideal + U2
    Vtheta2 = sigma * Vtheta2ideal
    VthetaRel2 = Vtheta2 - U2
    Vrel1 = (Vr1 ** 2 + VthetaRel1 ** 2) ** 0.5
    Vrel2 = (Vr2 ** 2 + VthetaRel2 ** 2) ** 0.5
    alpha1 = Vtheta1/Vr1
    alpha2 = Vtheta2/Vr2
    alphaRel1 = np.arctan(VthetaRel1/Vr1)
    alphaRel2 = np.arctan(VthetaRel2/Vr2)
    beta1 = alphaRel1 + inletangle * np.pi/180
    Vrel1 = (Vr1**2 + Vtheta1**2)**0.5
    Vrel2 = Vr2/(np.cos(alphaRel2))
    V1 = (Vr1)
    V2 = (Vr2 ** 2 + Vtheta2 ** 2) ** 0.5
    return(omega, Vtheta1, Vtheta2, VthetaRel1, VthetaRel2, Vr1, Vr2, alpha1, alpha2, alphaRel1, alphaRel2, beta1, beta2, V1, V2, Vrel1, Vrel2)

def bladenumber():
    # calculates minimum number of blades
    rmid = (r2 + r1) * 0.5
    betamid = np.arctan(0.5*(np.tan(beta1) + np.tan(beta2)))
    drVthetadr = (r2 * Vtheta2 - r1 * Vtheta1)/(r2 - r1)
    Vr = mdot/(2 * np.pi * density * rmid * h)
    Wavg = Vr/np.cos(betamid)
    Nb = mdot * drVthetadr / (2 * density * Wavg ** 2 * rmid * h)
    return(Nb)

def throatspace(r, theta):
    # finds minimum throat area in (ie min distance from any point to any point)
    mindist = 1000
    for i in range(0, len(r)):
        for j in range(0, len(r)):
            x1 = r[i] * np.cos(theta[0][i])
            y1 = r[i] * np.sin(theta[0][i])
            x2 = r[j] * np.cos(theta[1][j])
            y2 = r[j] * np.sin(theta[1][j])
            dist = ((x2 - x1) ** 2 + (y2 - y1) ** 2) **  0.5
            if dist < mindist:
                mindist = dist
    minarea = mindist * h
    return minarea

def maxthroat(r, theta):
    # finds throat area maximum (ie min distance of any point to trailing edge)
    mindist = 1000
    for j in range(0, len(r)):
        x1 = r[-1] * np.cos(theta[0][-1])
        y1 = r[-1] * np.sin(theta[0][-1])
        x2 = r[j] * np.cos(theta[1][j])
        y2 = r[j] * np.sin(theta[1][j])
        x3 = r[j] * np.cos(theta[-1][j]) # Extra point accounts for both forward sweep and backwards sweep
        y3 = r[j] * np.sin(theta[-1][j])
        dist1 = ((x2 - x1) ** 2 + (y2 - y1) ** 2) **  0.5
        dist2 = ((x3 - x1) ** 2 + (y3 - y1) ** 2) **  0.5
        if dist1 < mindist:
            mindist = dist1
        if dist2 < mindist:
            mindist = dist2
    minarea = mindist * h
    return minarea

        
def diffuserbladeplot(N, r1, r2, step, beta1, beta2):
    # plots diffuser blade shapes
    r = np.arange(r1, r2, 0.001)
    theta = [[0 for j in range(len(r))] for i in range(N)]
    fig, ax = plt.subplots(subplot_kw={'projection': 'polar'})
    beta = [0 for i in range(len(r))]
    for j in range(len(r)):
        beta[j] = beta1 + (beta2 - beta1) * (r[j] - r1)/(r2 - r1)
    for i in range(0, N): 
        theta[i][0] = i * 2 * np.pi/N
        for j in range(1, len(r)):
            theta[i][j] = theta[i][j-1] + np.tan(beta[j])/r[j] * step
        ax.plot(theta[i], r)
    ax.grid(True)
    plt.show()
    return(r, theta)
'''
betas = np.arange(-1.5, 1.5, 0.1)
omegas = []
velratios = []
for i in range(0, len(betas)):
    beta2 = betas[i]
    Vr2, omega, Vtheta1, Vtheta2, Vr1, alphaRel1, beta1, Vrel1, Vrel2 = velocitycalc()
    N = bladenumber()
    N = round(N * 1.25)
    Vr2, omega, Vtheta1, Vtheta2, Vr1, alphaRel1, beta1, Vrel1, Vrel2 = velocitycalc(N) # recalculates values using a refined blade number
    omegas.append(omega)
    velratios.append(Vrel2/Vrel1)
print(velratios)

'''
omega, Vtheta1, Vtheta2, VthetaRel1, VthetaRel2, Vr1, Vr2, alpha1, alpha2, alphaRel1, alphaRel2, beta1, beta2, V1, V2, Vrel1, Vrel2 = velocitycalc()
N = bladenumber()
N = round(N * 1.25)
omega, Vtheta1, Vtheta2, VthetaRel1, VthetaRel2, Vr1, Vr2, alpha1, alpha2, alphaRel1, alphaRel2, beta1, beta2, V1, V2, Vrel1, Vrel2 = velocitycalc(N) # recalculates values using a refined blade number
bladeplot(N, r1, r2, step, beta1, beta2)
print(omega, Vtheta1, Vtheta2, VthetaRel1, VthetaRel2, Vr1, Vr2, alpha1, alpha2, alphaRel1, alphaRel2, beta1, beta2, V1, V2, Vrel1, Vrel2)

r3 = r2 * vaneless
Vr3 = Vr2 * r2/r3
Vtheta3 = Vtheta2 * r2/r3
betad1 = np.arctan(Vtheta3/Vr3)
betad2 = 1.1

r, theta = diffuserbladeplot(10, r3, r4, step, betad1, betad2)
minarea = throatspace(r, theta)
exitarea = maxthroat(r, theta)
print("Area ratio is", exitarea/minarea)
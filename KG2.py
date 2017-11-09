"""
Eric Gelphman
University of California San Diego(UCSD)
Department of Mathematics
Irwin and Joan Jacobs School of Engineering Department of Electrical and Computer Engineering(ECE)
October 15, 2017

KrazyGlue script to numerically solve the mean curvature equation, a nonlinear second order ODE
Version 1.0.0
"""

import math
import matplotlib.pyplot as plt
from scipy.special import lambertw

def evalH(z, f, fp, fpp):
    r = 2*(lambertw((z**2 - f**2) / math.exp(1)).real) + 2.0
    H_Num = math.exp(r/4)*((r - 6)*(1 - fp**2)*(math.exp(-r/2)/(r - 1))*(fp*z - f) + r*fpp)
    H_Denom = 4*math.sqrt(2)*math.sqrt(r)*(1 - fp**2)**(3/2)
    return H_Num / H_Denom

def graphH(z_vals, H_vals):
    plt.plot(z_vals, H_vals)
    plt.xlabel("z")
    plt.ylabel("Mean Curvature H(z)")
    plt.show()

def main():
    z = 0.0
    f = 0.0
    fp = 0.0
    fpp = 0.0
    z_vals = []
    H_vals = []
    while(z < 10.0):
        if z > 0.0:
            f = 0.1 * math.sin(z)
        H = evalH(z, f, fp, fpp)
        z_vals.append(z)
        H_vals.append(H)
        z += 0.01
    graphH(z_vals, H_vals)


if __name__ == "__main__":
    main()

"""
Eric Gelphman
University of California San Diego(UCSD)
Department of Mathematics
Irwin and Joan Jacobs School of Engineering Department of Electrical and Computer Engineering(ECE)
September 17, 2017
KrazyGlue script to numerically solve the mean curvature equation, a nonlinear second order ODE
"""

import matplotlib.pyplot as plt
import fGenerator

#Function to graph r vs. H
def graphH(r_vals, H_vals):
    plt.plot(r_vals, H_vals)
    plt.xlabel("r")
    plt.ylabel("H(r)")
    plt.show()

def main():
    data = fGenerator.generator4()
    ofilen = open("output7.txt", 'r+')
    #Lists to store data to graph H
    r_vals = []
    H_vals = []
    H = 0.0
    #Traversal to iterate through data containing necessary information stored in 4-tuples
    #Step size is 1 index in list
    for i in data:
        r = i[0]
        fbar = i[1]
        dfbar_1 = i[2]
        dfbar_2 = i[3]
        if(r <= -40.0):
            H = fGenerator.fulmine.calcHLNr(r, fbar, dfbar_1, dfbar_2)#Pass values of r, fbar, fbar', fbar''
        else:
            H = fGenerator.fulmine.calcH(r, fbar + r - 4.0/3.0, dfbar_1 + 1, dfbar_2)#Pass values of r, f, f', f'' = fbar''
        line = "r: " + str(r) + " H: " + str(H) + "\n"
        ofilen.write(line)
        r_vals.append(r)
        H_vals.append(H)
    ofilen.close()
    graphH(r_vals, H_vals)

if __name__ == "__main__":
    main()

"""
Eric Gelphman
University of California San Diego(UCSD)
Department of Mathematics
Irwin and Joan Jacobs School of Engineering Department of Electrical and Computer Engineering(ECE)
January 16, 2018

KrazyGlue script to numerically solve the mean curvature equation, a nonlinear second order ODE
Version 1.1.0
"""

import matplotlib.pyplot as plt

#Function to read numerical data from file and to build array of 4-tuples
#(z, f(z) = w, f'(z), f''(z))
def readFromFile(fileName):
    file1 = open(fileName, 'r')
    list = []#Numerical values returned by fGenerator or calculateH, it works for both
    for line in file1:
        element = line.split(' ')
        list.append(element)
    return list

#Function to graph f(z)
def graphF(list):
    z_vals = []
    w_vals = []
    for tuple in list:
        z_vals.append(tuple[0])
        w_vals.append(tuple[1])
    plt.plot(z_vals, w_vals)
    plt.xlabel("z")
    plt.ylabel("f(z) = w")
    plt.show()

def graphH(list):
    z_vals = []
    H_vals = []
    for tuple in list:
        z_vals.append(tuple[0])
        H_vals.append(tuple[1])
    plt.plot(z_vals, H_vals)
    plt.xlabel("z")
    plt.ylabel("Mean Curvature H")
    plt.show()

def main():
    fList = readFromFile('fEvolve.txt')
    graphF(fList)
    #hList = readFromFile('outputH.txt')
    #graphH(hList)


if __name__ == "__main__":
    main()

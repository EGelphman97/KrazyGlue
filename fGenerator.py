"""
Eric Gelphman
University of California San Diego(UCSD)
Department of Mathematics
Jacobs School of Engineering Department of Electrical and Computer Engineering(ECE)
September 22, 2017

fGenerator script to generate values of r, fbar, fbar', and fbar'' and store them in such a
way that KrazyGlue can then efficiently calculate the mean curvature H
Version 1.0.7
"""

import math
import fulmine
import matplotlib.pyplot as plt

#Function to generate values for r, fbar = f - r + 4/3, fbar', and fbar'' in all 3 regions:
#parametric, gluing, and linear
#Note: All arithmetic operations are done by Fulmine
#Returns list of 4-tuples with values organized in this format [r, f(r), f'(r), f''(r)]
def generator4():
    r = -5.0
    f_bar = 0.0
    x = 1.0e-307
    #data is a list of 4-tuples that holds the 4 essential values: r, fbar = f - r + 4/3, fbar', and fbar''
    data = []
    #Loop starts from left side (r -> negative infinity) in the parametric region
    while r < 2.0:
        if r <= 0.1:
            rfbar_para = fulmine.peFbar(x)#Evaluate r and fbar at this value for x
            r = rfbar_para[0]#rfbar_para[1] = fbar
            fbar = rfbar_para[1]
            deriv_para = fulmine.peCalcDeriv(x, r)#Calculate fbar' and fbar'' for this value of x
            data.append((r, rfbar_para[1], deriv_para[0], deriv_para[1])) #Append 4-tuple to list
            #r is defined in terms of x, so increase x to increase r and continue loop
            if r < -200.0:
                x *= 1.15
            elif r > -200.0 and r < -5.0:
                x *= 1.03
            else:
                x += 1.0e-4
        else:
            fbar = fulmine.eeFbar(r)#Evaluate fbar explicitly in terms of r
            deriv_explicit = fulmine.eeCalcDeriv(r)#Calculate fbar' and fbar'' for this value of r
            data.append((r, fbar, deriv_explicit[0], deriv_explicit[1]))#Append 4-tuple to list
            r += 0.001#fbar is defined explicitly in terms of r here
    return data  #return list


#Function to write the numerical data to a file
def writeToFile(filename, data):
    file1 = open(filename, 'r+')#open file
    for i in data:#Write each 4-tuple to file
        file1.write("r: " + str(i[0]) + " fbar: " + str(i[1]) + " fbar': " + str(i[2]) + " fbar'': " + str(i[3]) + "\n")
    file1.close()#close file

#Function to graph fbar(r)
def graphf(r_vals, t_vals):
    plt.plot(r_vals, t_vals)
    plt.xlabel("r")
    plt.ylabel("fbar(r)")
    plt.show()

def main():
    fData = generator4()#Do entire numerical evaluation
    print(str(len(fData)))
    writeToFile("test7.txt", fData)#Write numerical data to file
    r_vals = []
    fbar_vals = []
    for i in fData:
        r_vals.append(i[0])
        fbar_vals.append(i[1])
    graphf(r_vals, fbar_vals)#Graph r vs. t

if __name__ == "__main__":
    main()
   

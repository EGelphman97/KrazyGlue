#Eric Gelphman
#June 27, 2017
#KrazyGlue Algorithm to investigate the mean curvature of spacetime near Black Holes
#fGenerator script to generate values of f
import math
import Fulmine
import matplotlib.pyplot as plt

#Function to generate values for r, fbar = f - r - 4/3, fbar', and fbar'' in all 3 regions:
#parametric, gluing, and linear
#Note: All arithmetic operations are done by Fulmine
#Returns list of 4-tuples with values organized in this format [r, f(r), f'(r), f''(r)]
def generator4():
    r = -5.0
    f_bar = fp = fpp = 0.0
    x = 1.0e-308
    #data is a list of 4-tuples that holds the 4 essential values: r, fbar = f - r + 4/3, fbar', and fbar''
    data = []
    #Traversal starts from left side (r -> negative infinity) in the parametric region
    while r <= -1.0:
        rfbar_para = peFbar(x)
        deriv_para = peCalcDeriv(x)
        data.append((rfbar_para[0], rfbar_para[1], deriv_para[0], deriv_para[1])) #Append to list
        print("x: " + str(x) + " r: " + str(r) + " f_bar: " + str(f_bar) + "\n")
        x += 0.000000000000001
    while r > -1.0 and r <= 2.0:
        fbar = eeFbar(r)
        deriv_explicit = eeCalcDeriv(r)
        data.append((r, fbar, deriv_explicit[0], deriv_explicit[1]))
        r += 0.01
    return data  #return list



def writeToFile(filename, data):
    file1 = open(filename, 'r+') #open file
    for i in data:
         #write to file
        file1.write("r: " + str(i[0]) + " f_bar: " + str(i[1]) + " f_bar': " + str(i[2]) + " f_bar'': " + str(i[3]) + "\n")
    file1.close()#close file

#Function to graph fbar(r)
def graphf(r_vals, t_vals):
    plt.plot(r_vals, t_vals)
    plt.xlabel("r")
    plt.ylabel("fbar(r)")
    plt.show()

def main():
    fData = generator4()
    #Ask user if they want file of numerical data in addition to graph
    print("Enter 1 to generate file of numerical data, enter 0 to end\n")
    response = int(input())
    if response == 1:
        writeToFile("test7.txt", fData)
    #Graph r vs. t
    r_vals = []
    fbar_vals = []
    for i in fData:
        r_vals.append(i[0])
        fbar_vals.append(i[1])
    graphf(r_vals, fbar_vals)

if __name__ == "__main__":
    main()

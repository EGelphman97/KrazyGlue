#Eric Gelphman
#June 27, 2017
#KrazyGlue Algorithm to investigate the mean curvature of spacetime near Black Holes
#fGenerator script to generate values of f
import math
import matplotlib.pyplot as plt

#Function to generate values for r, fbar = f - r - 4/3, fbar', and fbar'' in the parametric region
#r_init, f_bar_init, and fp_init are initial values for r, fbar, fbar', and fbar'' respectively
#Returns list of 4-tuples with values organized in this format [r, f(r), f'(r), f''(r)]
def parametricGenerator():
    r = -10002.0
    f_bar = fp = fpp = 0.0
    x = 0.5
    idx = 0
    #data is a list of 4-tuples that holds the 4 essential values: r, fbar = f - r + 4/3, fbar', and fbar''
    data = []
    #Traversal starts from left side (r -> negative infinity)
    while r <= -2.0:
        a = 1 - x
        r = 4/(3*a) + 4*a - 4*math.arctanh(a)
        #will need to do taylor series expansion of fbar for large negative values of r
        if r < -60.0:
            f_bar = -2*x + (5/2)*x**2 - (35/12)*x**3
        f_bar = 4/3 - 4/(3*a**(3/2))            #fbar = f - r + 4/3
        data.append((r, f_bar, 0.0, 0.0)) #Append to list
        if idx >= 2:
            #Assign values to index idx - 1 in list
            #Finite difference differentiation
            val_forward = data[idx+1]     #tuple at index idx + 1
            val_backward = data[idx-1]    #tuple at index 
            fp = (val_forward[1] - val_backward[1]) / (val_forward[0] - val_backward[0])   #fbar'
            fpp = (val_forward[1] - 2*f_bar + val_backward[1]) / (val_forward[0] - val_backward[0])**2        #fbar''
            data[idx-1][2] = fp
            data[idx-1][3] = fpp
        x += 0.01
        idx += 1    
    return data  #return list

#Function to generate values of r, fbar, fbar', and fbar'' for the glue and constant regions
#gr_lower is lower bound of glue region, gr_upper is upper bound for glue region, r_upper is highest r value in traversal
#Returns list of 4-tuples with values organized in this format [r, f(r), f'(r), f''(r)]
def otherGen(gr_lower, gr_upper, r_upper):
    data = []
    return data
    
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
    fData = parametricGenerator()
    fData.extend(otherGen(-1.0, -0.8, 2.0))
    #Graph r vs. t
    r_vals = []
    fbar_vals = []
    for i in fData:
        r_vals.append(i[0])
        fbar_vals.append(i[1])
    graphf(r_vals, fbar_vals)
    #Ask user if they want file of numerical data in addition to graph
    print("Enter 1 to generate file of numerical data, enter 0 to end\n")
    response = input()
    if response == 1:
        writeToFile("test5.txt", fData)


if __name__ == "__main__":
    main()

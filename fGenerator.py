#Eric Gelphman
#June 27, 2017 
#KrazyGlue Algorithm to investigate the mean curvature of spacetime near Black Holes
#fGenerator script to generate values of f
import math 
import matplotlib.pyplot as plt

#Function to generate values for r, f(r), f'(r), and f''(r)
#Filename is name of file where data will be stored
#Returns list of 4-tuples with values organized in this format [r, f(r), f'(r), f''(r)]
def parametricGenerator(filename):
    r = -10000.0
    t = 0.0
    y = 2520.0
    #data is a list of 4-tuples that holds the 4 essential values: r, t = f(r), f', and f''
    data = []
    file1 = open(filename, 'r+') #open file
    #Traversal starts from left side (r -> negative infinity)
    while r <= 2.0:
        #f(r) defined parametrically for r <= -1
        #where r = -1.0, y = 1.5451
        if r <= -1.20:
            a = math.tanh(y)
            t = 4*a - 4*y
            r = (4.0 / (3.0*(a**(3/2)))) + t
            #For large values of y, cosh(y) -> infinity, so sech(y)-> 0 (b terms) for large y
            #Case for when cosh(y) is not large
            try:
                b = 1 / ((math.cosh(y))**2)
                fp = (2*a**(9/2)) / (b + 2.0*a**(9/2))
                fpp = ((7.0 + 2.0*math.cosh(2.0*y))*(b**2)*(a**6)) / (((2.0*b + 4.0*a**(9/2)))**3)
            #Case for when cosh(y) is large    
            except OverflowError:
                fp = 1.0
                fpp = 0.0
            y -= 0.01
        #glue region
        elif r > -1.20 and r < -0.80:
            #g(r) = -0.8494r^3 - 3.6751r^2 - 4.2494r - 3.9823
            t = -0.8494*r**3 - 3.6751*r**2 - 4.2494*r - 3.9823
            fp = -2.5481*r**2 - 7.3503*r - 4.2494
            fpp = -5.0963*r - 7.3503
            r += 0.01  
        #Region where f(r) is constant    
        else:
            t = -2.50
            fp = 0.0
            fpp = 0.0
            r += 0.01
        data.append((r, t, fp, fpp))    #Append to list
        #write to file
        file1.write("r: " + str(r) + " t: " + str(t) + " f': " + str(fp) + " f'' " + str(fpp) + "\n")         
    file1.close() #close file
    return data  #return list   
     
#Function to graph f(r)    
def graphf(r_vals, t_vals):
    plt.plot(r_vals, t_vals)
    plt.xlabel("r")
    plt.ylabel("t")
    plt.show()    

def main():    
    filen = "test5.txt"
    data = parametricGenerator(filen)
    #Graph r vs. t
    i = 0
    r_vals = []
    t_vals = []
    while i < len(data):
        r_vals.append(data[i][0])
        t_vals.append(data[i][1])
        i += 1
    graphf(r_vals, t_vals)
    

if __name__ == "__main__":
    main()        
   


     
            
            
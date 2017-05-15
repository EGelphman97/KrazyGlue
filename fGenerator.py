#Eric Gelphman
#May 15, 2017 
#KrazyGlue Algorithm to investigate the mean curvature of spacetime near Black Holes
#fGenerator script to generate values of f
import math 
import matplotlib.pyplot as plt

def parametricGenerator(filename):
    s = 0.0 
    r = 0.0
    t = 0.0
    r_vals = []
    t_vals = []
    S_MAX = 10 - 2*(2/9)**(1/3)
    file1 = open(filename, 'r+')
    #Starting traversal from left side
    while s <= S_MAX:
        rho = (10 - s) * (9/2)**(1/3)
        t = 2*math.sqrt(2*rho) - 4*math.atanh(math.sqrt(2/rho))
        r = (10 - s)**(3/2) + t
        r_vals.append(r)
        t_vals.append(t)
        file1.write(str(r) + "_" + str(t) + "\n")
        s += 0.01
    #glue region
    while r < -3.0:
        #g(r) = -141.949r^3 - 1338.955r^2 - 4201.114r - 4388.363 
        t = -141.949*r**3 - 1338.955*r**2 - 4201.114*r - 4388.363
        r += 0.01
        r_vals.append(r)
        t_vals.append(t)
        file1.write(str(r) + "_" + str(t) + "\n")
    while r <= 2.0:
        file1.write(str(r) + "_-3.0\n")
        r_vals.append(r)
        t_vals.append(3.0)
        r += 0.01
    file1.close()
    return (r_vals, t_vals)     
        
def graphf(r_vals, t_vals):
    plt.plot(r_vals, t_vals)
    plt.xlabel("r")
    plt.ylabel("t")
    plt.show()    

def main():        
    #################################Start of Actual Script######################## 
    filen = "test4.txt"
    data = parametricGenerator(filen)
    graphf(data[0], data[1])

if __name__ == "__main__":
    main()        
   


     
            
            
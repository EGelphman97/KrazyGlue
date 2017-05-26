#Eric Gelphman
#May 26, 2017 
#KrazyGlue Algorithm to investigate the mean curvature of spacetime near Black Holes
#fGenerator script to generate values of f
import math 
import matplotlib.pyplot as plt

def parametricGenerator(filename):
    r = -10000.0
    t = 0.0
    s = 0.0 
    #data is a 4-tuple that holds the 4 essential values: r, t = f(r), f', and f''
    data = []
    file1 = open(filename, 'r+')
    while r <= 2.0:
        #Starting traversal from left side
        if r <= -3.2798:
            rho = (10 - s) * (9/2)**(1/3)
            t = 2*math.sqrt(2*rho) - 4*math.atanh(math.sqrt(2/rho))
            r = (10 - s)**(3/2) + t
            rho3 = rho**3 - 1
            fp = -5.45137/((2*rho)**(3/2)) - 3.85468/(rho**(3/2)*rho3) - 23.1281*rho**(3/2)/(rho3**2) + 3/(4*(10-s)**(1/2))
            fpp = -27.0001/((2*rho)**(5/2)) - 9.54588/(rho**(5/2)*rho3) + 38.1835*(rho**(1/2))/(rho3**2) - 229.101*rho**(7/2)/(rho3**3) + 3/(8*(10-s)**(3/2))
            data.append((r, t, fp, fpp))
            s += 0.001 
        #glue region
        elif r > -3.2798 and r < -3.0:
            #g(r) = -18.876r^3 - 179.157r^2 - 565.283r - 597.495
            t = -18.876*r**3 - 179.157*r**2 - 565.283*r - 597.495
            fp = -56.629*r**2 - 358.314*r - 565.283
            fpp = -113.257*r -358.314
            data.append((r, t, fp, fpp))
            r += 0.01  
        else:
            t = -4.4
            fp = 0.0
            fpp = 0.0
            data.append((r, t, fp, fpp))
            r += 0.01
        file1.write("r: " + str(r) + " t: " + str(t) + " f': " + str(fp) + " f'' " + str(fpp) + "\n")         
    file1.close()
    return data     
        
def graphf(r_vals, t_vals):
    plt.plot(r_vals, t_vals)
    plt.xlabel("r")
    plt.ylabel("t")
    plt.show()    

def main():        
    #################################Start of Actual Script######################## 
    filen = "test4.txt"
    data = parametricGenerator(filen)
    #graphf(data[0], data[1])

if __name__ == "__main__":
    main()        
   


     
            
            
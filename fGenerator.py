#Eric Gelphman
#May 15, 2017 
#KrazyGlue Algorithm to investigate the mean curvature of spacetime near Black Holes
#fGenerator script to generate values of f
import math 
import matplotlib.pyplot as plt

def parametricGenerator(filename):
    r = -10000.0
    t = 0.0
    r_vals = []
    t_vals = []
    file1 = open(filename, 'r+')
    #Starting traversal from left side
    while r <= -3.50:
        file1.write(str(r) + "_-5.0\n")
        r_vals.append(r)
        t_vals.append(-5.0)
        r += 0.01
    #glue region
    while r < -3.2798:
        #g(r) = -38.781r^3 - 392.681r^2 - 1323.548r - 1489.831
        t = -38.781*r**3 - 392.681*r**2 - 1323.548*r - 1489.831
        r_vals.append(r)
        t_vals.append(t)
        file1.write(str(r) + "_" + str(t) + "\n")  
        r += 0.01  
    S_START = 8.728016
    s = S_START    
    while r <= 2.0:
        rho = (10 - s) * (9/2)**(1/3)
        t = 2*math.sqrt(2*rho) - 4*math.atanh(math.sqrt(2/rho))
        r = (10 - s)**(3/2) + t
        r_vals.append(r)
        t_vals.append(t)
        file1.write(str(r) + "_" + str(t) + "\n")
        s -= 0.001  
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
   


     
            
            
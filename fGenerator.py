#Eric Gelphman
#May 14, 2017 
#KrazyGlue Algorithm to investigate the mean curvature of spacetime near Black Holes
#fGenerator script to generate values of f
import math 

def parametricGenerator(filename):
    s = 0.0 
    S_MAX = 10 - 2*(2/9)**(1/3)
    file1 = open(filename, 'r+')
    while s <= S_MAX:
        rho = (10 - s) * (9/2)**(1/3)
        t = 2*math.sqrt(2*rho) - 4*math.atanh(math.sqrt(2/rho))
        r = (10 - s)**(3/2) + t
        file1.write(str(r) + "_" + str(t) + "\n")
        s += 0.01
    file1.close()     
        
    

def main():        
    #################################Start of Actual Script######################## 
    filen = "test4.txt"
    parametricGenerator(filen)

if __name__ == "__main__":
    main()        
   


     
            
            
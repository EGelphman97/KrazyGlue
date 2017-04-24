#Eric Gelphman
#April 23, 2017 
#KrazyGlue Algorithm to investigate the mean curvature of spacetime near Black Holes

import matplotlib.pyplot as plt

#Perform a linear iterpolation if value of r is not found in numerical data
def linInterpolate(r, data, low, high):
    m = (data[high][1] - data[low][1]) / (data[high][0] - data[low][0]) #caluclate slope
    return m * (r - low[0]) + low[1]   #point-slope form
                
#Function to perform a binary search through a list of tuples representing 
#values of r and f(r)
#Note: In tuple structure, tuple[0] = r, tuple[1] = f(r)
def binarySearch(r, step_size, data, low, high):
    EPSILON = step_size / 10.0
    if low <= high:
        mid = int((low + high) / 2)
        if abs(data[mid][0] - r) <= EPSILON:         #r value is withtin epsilon
            return data[mid][1]
        else:
            if data[mid][0] - r < 0:            #data[mid][0] < r 
                return binarySearch(r, step_size, data, mid + 1, high)
            else:                                     #data[mid][0] > r
                return binarySearch(r, step_size, data, low, mid - 1)
    else:
        print("interpolation needed")
        return linInterpolate(r, data, low, high)

#Function to calculate H using necessary parameters       
def calculateH(r, t, fp, fpp): 
    X = 0.0 
    Y = 0.0
    dX_dr = 0.0
    dX_dt = 0.0  
    dY_dr = 0.0
    dY_dt = 0.0
    #Calculate X, Y, and their partial derivatives
    #Note: Need to account for where divide by 0 errors occur, as well as
    #when a negative value is under the square root 
    if r <= 0.9:                       #Schwarzchild Space
        X = 2 / ((6 * (r - t)) ** (1/3))
        Y = 1.65096 * ((r - t) ** (2/3)) 
        dX_dr = -0.36688 / ((r - t) ** (4/3))
        dX_dt = -dX_dr
        dY_dr = 1.10064 / ((r - t) ** (1/3))
        dY_dt = -dY_dr  
    elif r >= 1.1:                    #Friedman Space
        X = 1.65096 * ((1 - t) ** (2/3))
        Y = 1.65096 * (((r ** 3) * ((1 - t) ** 2)) ** (1/3))
        #dX_dr = 0.0
        dX_dt = (-1.10064 * (r ** 2)) / ((-((r ** 6)*(t - 1))) ** (1/3))
        dY_dr = (1.65096 * (r ** 2) * ((1 - t) ** 2)) / (((r ** 3) * ((1 - t) ** 2)) ** (2/3))
        dY_dt = (-1.10064 * (r ** 3) * (1 - t)) / (((r ** 3) * ((1 - t) ** 2)) ** (2/3))
    else:
        X_num = (2.4825 * r) - (1.655 * t) + 0.2553
        X_denom = (8.2171 * r ** 3) + (2.5346 * r ** 2) - (16.4342 * r ** 2 * t) - (3.6559 * r) + (9.7215 * r * t) - 1.4377 * t + 0.6469
        X = X_num / (X_denom ** (1/3))
        Y = ((3.7238 * r ** 2) - (7.4475 * r * t) + 2.25 * r + 2.2028 * t - 0.9912) ** (1/3) 
        dX_dr = (2.4825 / (X_denom ** (1/3))) - ((X_num * (24.6513 * r ** 2) - (32.8684 * r * t) - 5.0692 * r -3.6559) / (3 * (X_denom ** (4/3))))
        dX_dt = (((1.4377 + 16.4342 * r ** 2) * X_num) / (3 * (X_denom ** (4/3)))) - (1.655 / (X_denom ** (1/3)))
        dY_dr = (2.25 + 7.74476 * r) / (3 * (((3.7238 * r ** 2) - (7.4475 * r * t) + 2.25 * r + 2.2028 * t - 0.9912) ** (2/3)))
        dY_dt = 0.7343 / (((3.7238 * r ** 2) - (7.4475 * r * t) + 2.25 * r + 2.2028 * t - 0.9912) ** (2/3))
    #Check so a negative number is not under the square root in the equation for H
    if abs(X) < abs(fp):
        return "not real"
    else:    
        H_Num = (2*fp*dY_dt) - ((X**3)*Y*dX_dr) + (X*Y*fp*(dX_dt + 2*fp*dX_dr)) - (2*(X**4)*dY_dr) - ((X**2)*(Y*fpp+2*fp*(dY_dt -fp*dY_dr)))
        H_Denom = Y * ((X**2)/((X**2)-(fp**2))**(1/2)) * (((X**2)- (fp**2))**2)
    #Check whether H is undefined
    if H_Denom == 0.0:
        return "undefined"
    else:
        return H_Num / H_Denom
    
#Function to read values of r, f(r) from file and store them as a tuple <r, f(r)>
#The tuples are then stored in ascending r-order in a list, which is returned   
def readFromFile():
    with open("test3.txt", "r") as lines:
        f_data = []
        for line in lines:
            split = line.partition("_")
            try:
                f_data.append(((round(float(split[0]), 4)), round(float(split[2]), 5)))
            except ValueError:
                break
    lines.close()        
    return f_data
    
#Function to check boundary conditions around r = 1.5   
#f'(r) and f''(r) need to be 0 at r = 1.5 or about there 
def checkBoundaryCondition(r, f, fp, fpp, EPSILON):
    if r >= 1.4 and r <= 1.6:
        if(abs(f - r + 4/3) >= EPSILON or abs(fp - 1) >= EPSILON  or abs(fpp) >= EPSILON):
            print("Boundary conditions not met!")
            print("r = " + str(r) + " f(r) = " + str(f) + " f'(r) = " + str(fp) + " f''(r) = " + str(fpp))    
        else:
            print("Boundary Conditions Met!")  

def graphH(r_vals, f_vals, H_vals):
    #plt.plot(r_vals, f_vals)
    #plt.xlabel("r")
    #plt.ylabel("f(r)")
    #plt.show()
    plt.plot(r_vals, H_vals)
    plt.xlabel("r")
    plt.ylabel("H(r)")
    plt.show()
    
    
            
            
        
def main():     
    STEP_SIZE = 0.1
    f_data = readFromFile()
    file2 = open("output3.txt", 'r+')
    r_vals = []
    f_vals = []
    H_vals = []
    #Traversal
    r = 2.0
    f = -10.0
    while r >= -10000:
         #f_data must be sorted in increasing order
        f = binarySearch(r, STEP_SIZE, f_data, 0, len(f_data) - 1)
        if r >= -9.9:      #Not in the glue region
            fp = 0.0
            fpp = 0.0
        elif r <= -10.1:
            fp = 1 + 1/(3*r**2) - 8/(9*r**3)
            fpp = -2/(3*r**3) + 8/(3*r**4)
        else:
            fp = -1007.7*r**2 - 20159.02*r - 100809.622
            fpp = -2015.4*r -20159.02
        #Check boundary conditions around r = 1.5
        #if (r >= 1.4 and r <= 1.6):
            #checkBoundaryCondition(r, f, fp, fpp, EPSILON)  
        output = calculateH(r, f, fp, fpp)
        if (output != "not real" and output != "undefined"):    
            H = output
            line = "r: " + str(r) + " f: " + str(f) + " fp: " + str(fp) + " fpp: " + str(fpp) + " H: " + str(H) + "\n"
            file2.write(line)
            r_vals.append(r)
            f_vals.append(f)
            H_vals.append(H)
        if r < -1:
            STEP_SIZE += 0.1
        r -= STEP_SIZE
    graphH(r_vals, f_vals, H_vals)        
        
        
if __name__ == "__main__":
    main()        
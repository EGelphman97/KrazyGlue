#Eric Gelphman
#April 17, 2017 
#KrazyGlue Algorithm to investigate the mean curvature of spacetime near Black Holes

#Perform a linear iterpolation if value of r is not found in numerical data
def linInterpolate(r, data, low, high):
    m = (data[high][1] - data[low][1]) / (data[high][0] - data[low][0]) #caluclate slope
    return m * (r - low[0]) + low[1]   #point-slope form
                
#Function to perform a binary search through a list of tuples representing 
#values of r and f(r)
#Note: In tuple structure, tuple[0] = r, tuple[1] = f(r)
def binarySearch(r, data, low, high):
    EPSILON = 0.01
    if low <= high:
        mid = int((low + high) / 2)
        if abs(data[mid][0] - r) <= EPSILON:         #r value is withtin epsilon
            return data[mid][1]
        else:
            if data[mid][0] - r < 0:            #data[mid][0] < r 
                return binarySearch(r, data, mid + 1, high)
            else:                                     #data[mid][0] > r
                return binarySearch(r, data, low, mid - 1)
    else:
        print(str(data[mid][0]) + " " + str(data[mid][1]))
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
        X = 2 / (6 * (r - t) ** (1/3))
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
        
def main():     
    STEP_SIZE = 0.1
    f_data = readFromFile()
    file2 = open("output3.txt", 'r+')
    #Traversal
    r = 2.0
    f_prev = -10.0
    f = f_prev
    fprime_prev = 0.0
    r_prev = 2.1
    while r >= -9999.7:
        #f_data must be sorted in increasing order
        f = binarySearch(r, f_data, 0, len(f_data) - 1)
        fp = f - f_prev / (r - r_prev)
        fpp = fp - fprime_prev / (r - r_prev)
        #Check boundary conditions around r = 1.5
        #if (r >= 1.4 and r <= 1.6):
            #checkBoundaryCondition(r, f, fp, fpp, EPSILON)  
        output = calculateH(r, f, fp, fpp)   
        print(output)
        if (output != "not real" and output != "undefined"):    
            H = output
            line = "r: " + str(r) + " f: " + str(f) + " fp: " + str(fp) + " fpp: " + str(fpp) + " H: " + str(H) + "\n"
            file2.write(line)
            f_prev = f
            fprime_prev = fp  
            r_prev = r
        if r < -1:
            STEP_SIZE += 0.1
        r -= STEP_SIZE
        
if __name__ == "__main__":
    main()        
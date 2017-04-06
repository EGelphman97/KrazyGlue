#Eric Gelphman
#February 9, 2017 
#KrazyGlue Algorithm to investigate the mean curvature of spacetime near Black Holes

#Function to perform a binary search through a list of tuples representing 
#values of r and f(r)
#Note: In tuple structure, tuple[0] = r, tuple[1] = f(r)
def binarySearch(r, data):
    low = 0
    high = len(data) - 1
    found = False
    while low <= high and not found:
        pos = 0
        mid = int((low + high) / 2)
        if abs(data[mid][0] - r) <= 0.05:
            pos = mid
            found = True
        else:
            if abs(data[mid][0] - r) > 0.05:
                high = mid - 1
            else:
                low = mid + 1
    if found:                      #Value of f found
        return data[pos][1]
    else:                          #Value of f not found
        if low == high:
            high += 1
        #linear interpolation
        m = (data[high][1] - data[low][1]) / (data[high][0] - data[low][0])            
        return (m * r) - (m * data[high][0]) + data[high][1]

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
    # TODO: Automtaic Gluing
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
def readFromFile(file1):
    f_data = []
    for line in file1:
        split = line.partition("_")
        f_data.append((float(split[0]), float(split[2])))
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
    STEP_SIZE = float(input("Enter the step size: "))
    file1 = open(input("Enter the name of the file that contains values of r and f(r): "), 'r')
    f_data = readFromFile(file1)
    file1.close()
    outputFileName = input("Enter the name of the file that will store the output values: ")
    file2 = open(outputFileName, 'r+')
    
    #Traversal
    r = 2.0
    f_prev = 0.0
    fprime_prev = 0.0
    r_prev = 0.0
    while r >= -10000:
        #f_data must be sorted in increasing order
        f = binarySearch(r, f_data)
        fp = f - f_prev / (r - r_prev)
        fpp = fp - fprime_prev / (r - r_prev)
        #Check boundary conditions around r = 1.5
        #if (r >= 1.4 and r <= 1.6):
            #checkBoundaryCondition(r, f, fp, fpp, EPSILON)
        output = calculateH(r, f, fp, fpp)    
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


    

        
    



    
        
    
    

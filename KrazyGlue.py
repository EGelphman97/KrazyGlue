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
    a = (4.5) ** (2/3)
    b = r - t 
    c = 1 - t
    pwr_13 = 1 / 3
    pwr_43 = 4 / 3
    #Calculate X, Y, and their partial based on r
    #Note: Need to account for where divide by 0 errors occur, as well as
    #when a negative value is under the square root 
    if r <= 0.9:                       #Schwarzchild Space
        X = 2 / ((6 * b) ** pwr_13)
        Y = (4.5 * (b ** 2)) ** pwr_13
        dX_dr = -4 / ((6 * b) ** pwr_43)
        dX_dt = 4 / ((6 * b) ** pwr_43)
        dY_dr = 3 / (a * b ** pwr_13)
        dY_dt = -3 / (a * b ** pwr_13)  
    elif r >= 1.1:                    #Friedman Space
        X = (3 * (r ** 2) * c) / (((6 * (r ** 6)) * c) ** pwr_13)
        Y = ((4.5 * r ** 3) * c ** 2) ** pwr_13
        dX_dr = (36 * (r ** 2 - (r ** 2 * t) - 1 + t)) / ((6 ** pwr_43) * (r ** 3) * (c ** pwr_13))
        dX_dt = (-12 * (r ** 8)) / (((6 * (r ** 6)) ** pwr_43) * c)
        dY_dr = (4.5 ** pwr_13) * (c ** (2/3)) 
        dY_dt = (-3 * r) / (a * (c ** pwr_13))   
    else:                             #Glue Region
        M = 8.0*r**3 - 14.925*r**2 + 7.425*r + 0.5747
        Mp = 24.0*r**2 - 29.85*r + 7.425
        t_0t = 0.5*r + 0.45 - t
        X = ((t_0t)*(M)**2)**(-1 * pwr_13)*(4.40256966519284*r**3 - 8.21354403162538*r**2 + 4.0861349705071*r + 0.550321208149104*(t_0t)*(Mp) + 0.31626959832329)
        Y = ((t_0t)**2*(4.5 * M))**pwr_13                             
        dX_dr = ((t_0t)*(M)**2)**(-1 * pwr_13)*((t_0t)*(19.8115634933678*r**2 - 24.6406320948762*r + 0.550321208149104*(48.0*r - 29.85)*(t_0t) + 6.12920245576065)*(M) - (pwr_43*r**3 - 2.4875*r**2 + 1.2375*r + pwr_13*(t_0t)*(2*M) + 0.0957833333333333)*(4.40256966519284*r**3 - 8.21354403162538*r**2 + 4.0861349705071*r + 0.550321208149104*(t_0t)*(Mp) + 0.31626959832329))/((t_0t)*(M))
        dX_dt = ((t_0t)*(M)**2)**(-1 * pwr_13)*(1.46752322173095*r**3 - 2.73784801054179*r**2 + 1.36204499016903*r + (t_0t)*(-13.2077089955785*r**2 + 16.4270880632508*r - 4.0861349705071) + 0.183440402716368*(t_0t)*(Mp) + 0.105423199441097)/(t_0t)
        dY_dr = ((t_0t)**2*(4.5 * M))**pwr_13*(pwr_13*(t_0t)*(4.5 * M) + (t_0t)**2*(36.0*r**2 - 44.775*r + 11.1375))/((t_0t)**2*(4.5 * M))
        dY_dt = ((t_0t)**2*(4.5 * M))**pwr_13*(-1 * pwr_13*r + (2/3)*t - 0.3)/(t_0t)**2
    H_Num = (2*fp*dY_dt) - ((X**3)*Y*dX_dr) + (X*Y*fp*(dX_dt + 2*fp*dX_dr)) - (2*(X**4)*dY_dr) - ((X**2)*(Y*fpp+2*fp*(dY_dt -fp*dY_dr)))
    H_Denom = Y * ((X**2)/((X**2)-(fp**2))**(1/2)) * (((X**2)- (fp**2))**2)
    #Check whether H is undefined
    try:
        H = H_Num / H_Denom
    except ZeroDivisionError:
        print("Error! H(r) undefined at r = ", r)
        quit()
    return H  

#Function to read values of r, f(r) from file and store them as a tuple <r, f(r)>
#The tuples are then stored in ascending r-order in a list, which is returned   
def readFromFile(file1):
    f_data = []
    for line in file1:
        split = line.partition("_")
        f_data.append((float(split[0]), float(split[2])))
    return f_data
    
#Function to check boundary conditions around r = 1.5   
def checkBoundaryCondition(r, f, fp, fpp, EPSILON):
    #f'(r) and f''(r) need to be 0 at r = 1.5 or about there  
    if(abs(f - r + 4/3) >= EPSILON or abs(fp - 1) >= EPSILON  or abs(fpp) >= EPSILON):
        print("Boundary conditions not met!")
        print("r = " + str(r) + " f(r) = " + str(f) + " f'(r) = " + str(fp) + " f''(r) = " + str(fpp))    
    else:
        print("Boundary Conditions Met!")    
        
    
def main():    
    #########################Actual Script Starts Here#############################    
    EPSILON = 0.00001
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
        if(r >= 1.4 and r <= 1.6):
            checkBoundaryCondition(r, f, fp, fpp, EPSILON)
        H = calculateH(r, f, fp, fpp)
        line = str(r) + "_" + str(f) + "_" + str(fp) + "_" + str(fpp) + "_" + str(H) + "\n"
        file2.write(line)
        f_prev = f
        fprime_prev = fp  
        r_prev = r
        if r < -1:
            STEP_SIZE += 0.1
        r -= STEP_SIZE
        
if __name__ == "__main__":
    main()        


    

        
    



    
        
    
    

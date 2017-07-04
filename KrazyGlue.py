#Eric Gelphman
#June 28, 2017 
#KrazyGlue Algorithm to investigate the mean curvature of spacetime near Black Holes

import matplotlib.pyplot as plt
import fGenerator

"""
#Perform a linear iterpolation if value of r is not found in numerical data
def linInterpolate(r, data, low, high):
    m = (data[high][1] - data[low][1]) / (data[high][0] - data[low][0]) #caluclate slope
    return m * (r - low[0]) + low[1]   #point-slope form
                
#Function to perform a binary search through a list of tuples representing 
#values of r and f(r)
#Note: In tuple structure, tuple[0] = r, tuple[1] = f(r)
def binarySearch(r, step_size, data, low, high):
    EPSILON = step_size / 10.0
    if EPSILON > 0.01:
        EPSILON = 0.01
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
        print("interpolation needed r: " + str(r) + "low: " + str(low) + "high: " + str(high))
        return linInterpolate(r, data, low, high)
"""        

#Function to calculate H given r, t = f(r), f'(r), and f''(r)
#Returns value of H for the 4-tuple
#Note: This value represents the solution to the 2nd Order Nonlinear ODE that
#represents H at the r-value given in the parameters  
def calculateH(r, t, fp, fpp): 
    X = 0.0 
    Y = 0.0
    dX_dr = 0.0
    dX_dt = 0.0  
    dY_dr = 0.0
    dY_dt = 0.0
    #Calculate X, Y, and their partial derivatives
    #Note: Need to account for where divide by 0 errors occur, as well as when a negative value is under the square root 
    #Schwarzchild Space: M(r) = 1, t0(r) = r, applies for all r <= 0.8
    if r <= 0.8:                       
        X = 2 / ((6 * (r - t)) ** (1/3))
        Y = 1.65096 * ((r - t) ** (2/3)) 
        dX_dr = -0.36688 / ((r - t) ** (4/3))
        dX_dt = -dX_dr
        dY_dr = 1.10064 / ((r - t) ** (1/3))
        dY_dt = -dY_dr 
    #Friedman Space: M(r) = r^3, t0(r) = 1, applies for all r >= 1.2   
    elif r >= 1.2:                    
        X = (3**0.6666666666666666*r**2*(1 - t)) / (2**0.3333333333333333*(r**6*(1 - t))**0.3333333333333333)
        Y = 1.6509636244473134*(r**3*(1 - t)**2)**0.3333333333333333
        dX_dr = -((6**0.6666666666666666*r**7*(1 - t)**2)/(r**6*(1 - t))**1.3333333333333333) + (6**0.6666666666666666*r*(1 - t))/(r**6*(1 - t))**0.3333333333333333
        dX_dt = (r**8*(1 - t))/(6**0.3333333333333333*(r**6*(1 - t))**1.3333333333333333) - (3**0.6666666666666666*r**2)/(2**0.3333333333333333*(r**6*(1 - t))**0.3333333333333333)
        dY_dr = (1.6509636244473134*r**2*(1 - t)**2)/(r**3*(1 - t)**2)**0.6666666666666666
        dY_dt = (-1.100642416298209*r**3*(1 - t))/(r**3*(1 - t)**2)**0.6666666666666666
    #Glue Region: M(r) = 4.25r^3 - 7.35r^2 + 3.6r + 0.648, t0(r) = -1.25r^2 + 3r - 0.8, applies for 0.8 < r < 1.2
    #Glue Region is for the gluing of the metric, and is different than gluing for f in fGenerator
    else:
        X = (2*(3 - 2.5*r)*(0.648 + 3.6*r - 7.35*r**2 + 4.25*r**3) + (3.6 - 14.7*r + 12.75*r**2)*(-0.8 + 3*r - 1.25*r**2 - t))/(6**0.3333333333333333*((0.648 + 3.6*r - 7.35*r**2 + 4.25*r**3)**2*(-0.8 + 3*r - 1.25*r**2 - t))**0.3333333333333333)
        Y = 1.6509636244473134*((0.648 + 3.6*r - 7.35*r**2 + 4.25*r**3)*(-0.8 + 3*r - 1.25*r**2 - t)**2)**0.3333333333333333	
        dX_dr = (0.3147944506877429 + 221.64186658205165*r**10 - 27.28675990405975*r**11 + r**9*(-768.8853224449925 - 28.066381615604318*t) + 0.573283980986305*t + 
        0.28691491788233897*t**2 + r**8*(1491.9154348970615 + 157.28180128901397*t) + 
        r**2*(-11.748431131370616 - 23.240019899778524*t - 6.166611143524113*t**2) + 
        r**4*(-12.129112454495935 + 56.34030249795842*t - 4.509139271939641*t**2) + 
        r**6*(1207.2353490473233 + 481.2279501109907*t - 1.2880880741155186e-14*t**2) + 
        r**7*(-1746.1824316337115 - 377.4219384095341*t - 1.6101100926443983e-15*t**2) + 
        r**5*(-407.8425333707173 - 314.37502635182875*t + 0.5609391702828131*t**2) + 
        r*(-0.1599506126960389 + 0.32576916425958746*t + 1.054371573642127*t**2) + 
        r**3*(53.00757425014704 + 47.767873798856755*t + 8.477385842364002*t**2))/((0.15247058823529414 + 0.8470588235294118*r - 1.7294117647058822*r**2 + 1.*r**3)**2*((0.648 + 3.6*r - 7.35*r**2 + 4.25*r**3)**2*(-0.8 + 3*r - 1.25*r**2 - t))**0.3333333333333333*(0.64 - 2.4*r + 1.*r**2 + 0.8*t))
        dX_dt = (-0.032915043135895446 + 12.877516270689048*r**9 - 1.55924342308913*r**10 + 
        r**6*(-73.09910285815452 - 33.51184234207549*t) + 
        r**4*(-16.593132364449772 - 16.986886983224093*t) + 
        r**8*(-41.33689197258338 - 3.742184215413911*t) + 
        r*(-0.27442848982705526 - 0.17262657366047984*t) + 
        r**2*(0.5558242093947334 + 0.8265530168350382*t) - 0.024563465026787648*t + 
        r**3*(2.5333600414330157 + 2.627242009102273*t) + 
        r**7*(71.25533105645985 + 17.258073087555918*t) + 
        r**5*(45.694958513673356 + 33.69090048686677*t))/((0.15247058823529414 + 0.8470588235294118*r - 1.7294117647058822*r**2 + 1.*r**3)**2*((0.648 + 3.6*r - 7.35*r**2 + 4.25*r**3)**2*(-0.8 + 3*r - 1.25*r**2 - t))**0.3333333333333333*(0.64 - 2.4*r + 1.*r**2 + 0.8*t))
        dY_dr = (-0.44377902225143756 - 143.16950180754043*r**5 + 25.581337410056026*r**6 + 
        r**3*(-280.0859788874867 - 96.58137203016783*t) + 1.030201301655124*t + 
        1.9811563493367759*t**2 + r**4*(295.79764938014364 + 29.23581418292117*t) + 
        r*(-16.351143736526197 - 34.93439029330516*t - 8.089721759791836*t**2) + 
        r**2*(120.09109404229754 + 98.89272110439406*t + 7.016595403901082*t**2))/((0.648 + 3.6*r - 7.35*r**2 + 4.25*r**3)*(0.8 - 3*r + 1.25*r**2 + t)**2)**0.6666666666666666
        dY_dt = (5.847162836584234*(0.15247058823529414 + 0.8470588235294118*r - 
        1.7294117647058822*r**2 + 1.*r**3)*(0.64 - 2.4*r + 1.*r**2 + 0.8*t))/((0.648 + 3.6*r - 7.35*r**2 + 4.25*r**3)*(0.8 - 3*r + 1.25*r**2 + t)**2)**0.6666666666666666
    #Check so a negative number is not under the square root in the equation for H
    if abs(X) < abs(fp):
        return "not real"
    else:    
        H_Num = (2*(fp**3)*dY_dr) - ((X**3)*Y*dX_dt) + ((X*Y*fp)*(dX_dr + 2*fp*dX_dt)) - (2*(X**4)*dY_dt) - ((X**2)*((Y*fpp)+(2*fp)*(dY_dr - fp*dY_dt)))
        H_Denom = X * Y * ((X**2 - fp**2)**(3/2))
    #Check whether H is undefined
    if H_Denom == 0.0:
        return "undefined"
    else:
        return H_Num / H_Denom


#Function to graph r vs. H            
def graphH(r_vals, H_vals):
    plt.plot(r_vals, H_vals)
    plt.xlabel("r")
    plt.ylabel("H(r)")
    plt.show() 
        
def main():    
    fileni = "test5.txt"
    data = fGenerator.parametricGeneratorparametricGenerator(fileni)
    fileno = open("output5.txt", 'r+')
    #Lists to store data to graph H
    r_vals = []
    H_vals = []
    #Traversal to iterate through data containing necessary information stored in 4-tuples
    #Step size is 1 index in list
    for i in data:
        r = data[i][0]
        f = data[i][1]
        fp = data[i][2]
        fpp = data[i][3]
        output = calculateH(r, f, fp, fp)
        if (output != "not real" and output != "undefined"):    
            H = output
            line = "r: " + str(r) + " f: " + str(f) + " fp: " + str(fp) + " fpp: " + str(fpp) + " H: " + str(H) + "\n"
            fileno.write(line)
            r_vals.append(r)
            H_vals.append(H)
    fileno.close()        
    graphH(r_vals, H_vals)        
        
        
if __name__ == "__main__":
    main()        
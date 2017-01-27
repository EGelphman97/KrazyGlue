#Eric Gelphman
#January 26, 2017 
#KrazyGlue Algorithm to investigate the mean curvature of spacetime near Black Holes
#Function class to handle evaluation and differentiation
#Note: Class has been modified to handle all operations for f as well as X and Y

from sympy import Symbol, Limit, S

class Function:
    #Constructor
    def _init_(self, expression1, h):
        self.expression = expression1
        self.r = Symbol('r')
        self.data = [[]]
        self.step_size = h
        self.EPSILON = 0.00001
    
    #Evaluate function f(r) at r = a    
    def evaluate(self, a):
        return self.expression.subs({self.r:a})
        
    #Numerically differentiate f(r) at r = a
    #Order parameter: order of derivative(first or second)
    def differentiate(self, a, order):
        h = 0.0000001
        if order == 1:
            return (self.expression.subs({self.r:a+h}) + self.expression.subs({self.r:a-h})) / (2 * h)
        elif order == 2:
            return (self.expression.subs({self.r:a+h}) - (2 * self.expression.subs({self.r:a})) + self.expression.subs({self.r:a-h})) / (h**2)
            
    #Generate a table of values for r, f(r), f'(r), f''(r), X, Y, and the partial derivatives of X and Y       
    def generateTable(self):
        if self.checkBoundaryCondition() == False :
            print("Boundary Condition 2(and Therefore 3) not met as f does not approach negative infinity")
            quit()
        r = 2.0
        t, X, Y, dX_dr, dX_dt, dY_dr, dY_dt = 0.0
        a = (4.5) ** (2/3)
        b = r - t 
        c = 1 - t
        pwr_13 = 1 / 3
        pwr_43 = 4 / 3
        while r >= -10000:
              t = self.evaluate(r)
              a = (4.5) ** (2/3)
              b = r - t 
              c = 1 - t
              pwr_13 = 1 / 3
              pwr_43 = 4 / 3
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
              self.data.append([r, t, self.differentiate(r, 1), self.differentiate(r, 2), X, Y, dX_dr, dX_dt, dY_dr, dY_dt])
              if r < -1:
                  self.step.size += 0.5
              r -= self.step_size
        if abs(self.differentiate(r, 2)) >= self.EPSILON:
            print("Boundary Condition 1(and Therefore 3) not met as f double prime does not equal 0")
            quit()
        return self.data  
    
    #Checks boundary condition that f(r) must go to negative infinity  
    def checkBoundaryCondition(self):
        lim = Limit(self.expression, self.r, S.Negative_Infinity)    
        if lim != S.NegativeInfinity: 
            print("Boundary Condition 2(and Therefore 3) not met as f does not approach negative infinity")
            return False
        else:
            return True
            
            
                
            
        

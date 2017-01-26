#Eric Gelphman
#January 16, 2017 
#KrazyGlue Algorithm to investigate the mean curvature of spacetime near Black Holes
#Function class to handle evaluation and differentiation
from sympy import Symbol, Limit, S

class Function_SVar:
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
            
    #Generate a table of values for r, f(r), f'(r), f''(r)        
    def generateTable(self):
        r = 2.0
        while r >= -10000:
            self.data.append([r, self.evaluate(r), self.differentiate(r, 1), self.differentiate(r, 2)])
            r -= self.step_size
        if abs(self.differentiate(r, 1)) >= self.EPSILON & abs(self.differentiate(r, 2)) >= self.EPSILON:
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
            
            
                
            
        

#Eric Gelphman
#January 16, 2017 
#KrazyGlue Algorithm to investigate the mean curvature of spacetime near Black Holes
#Function class to handle evaluation and differentiation
from sympy import Symbol

class Function:
    def _init_(self, expression1):
        self.expression = expression1
    def evaluate(self, a):
        r = Symbol('r')
        return self.expression.subs({r:a})
    def differentiate(self, a, order):
        r = Symbol('r')
        h = 0.0000001
        if order == 1:
            return (self.expression.subs({r:a+h}) + self.expression.subs({r:a-h})) / (2 * h)
        elif order == 2:
            return (self.expression.subs({r:a+h}) - (2 * self.expression.subs({r:a})) + self.expression.subs({r:a-h})) / (h**2)
            
                
            
        

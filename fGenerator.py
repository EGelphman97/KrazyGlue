#Eric Gelphman
#April 17, 2017 
#KrazyGlue Algorithm to investigate the mean curvature of spacetime near Black Holes
#fGenerator script to generate values of f
from sympy import Symbol, sympify, simplify

#Class definition to more easily handle evaluation and writing to file
class FunctionStruct:
    #Constructor
    def __init__(self, expression1, h, lbound, hbound):
        self.expression = expression1
        self.r = Symbol('r')
        self.step_size = h
        self.lower_limit = lbound
        self.upper_limit = hbound
    
    #Evaluate function f(r) at r = a    
    def generateFile(self, filename):
        r1 = self.lower_limit
        file1 = open(filename, 'r+')
        while r1 <= self.upper_limit:
            f = simplify(self.expression.subs({self.r:r1}))
            file1.write(str(r1) + "_" + str(f) + "\n")
            r1 += self.step_size
        file1.close()  


def main():        
    #################################Start of Actual Script########################        
    wfileN1 = "test3.txt"
    wfileN2 = "test3_1.txt"
    wfileN3 = "test3_2.txt"
    h = 0.1
    expr_f = sympify("r - 4/3 - 1/(3*r) + 4/(9*r**2)", evaluate=False)
    obj = FunctionStruct(expr_f, h, -10000, -10.2)
    obj.generateFile(wfileN1)
    expr_glue = sympify("342.850*r**3 + 10282.990*r**2 + 102975.017*r + 342491.827", evaluate=False)
    obj = FunctionStruct(expr_glue, h, -10.1, -9.9)
    obj.generateFile(wfileN2)
    r = -9.8
    fileW = open(wfileN3, 'r+')
    while r <= 2.0:
        fileW.write(str(r) + "_" + "-10.0" + "\n")
        r += h

if __name__ == "__main__":
    main()        
   


     
            
            
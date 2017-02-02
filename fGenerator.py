#Eric Gelphman
#February 1, 2017 
#KrazyGlue Algorithm to investigate the mean curvature of spacetime near Black Holes
#fGenerator script to generate values of f
from sympy import Symbol

#Class definition to more easily handle evaluation and writing to file
class FunctionStruct:
    #Constructor
    def _init_(self, expression1, h, lbound, hbound):
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
            f = self.expression.subs({self.r:r1})
            line = str(r1) + "_" + str(f)
            file1.write(line)
            r1 += self.step_size
        file1.close()  
        
#################################Start of Actual Script########################        
wfileN = input("Enter name of file: ")
expr = input("Enter and expression f(r): ")
h = input("Enter the step size: ")
lbound = input("Enter the lowest r value for which this segment is defined")
hbound = input("Enter the largest r value for which this segment is defined")
obj = FunctionStruct(expr, h, lbound, hbound)
obj.generateFile(wfileN)
   


     
            
            
/*Eric Gelphman
  University of California San Diego(UCSD)
  Department of Mathematics
  Jacobs School of Engineering Department of Electrical and Computer Engineering(ECE)
  September 4, 2017

  Fulmine, a C subroutine to perform numerical differentiation and evaluation
  for fGenerator and KrazyGlue.
*/

//#include <Python.h>//Library needed to interface with Python
#include <stdio.h>
#include <math.h>

/*Function to evaluate r and fbar in all 3 regions: parametric (r <= -1.0), glue region
(-1.0 < r < -0.6), and linear(r >= -0.6)
*/
double* rfEval(double x)
{

    double a,r,fbar;
    double result[2];
    double* resultPtr;
    resultPtr = result;
    //Evaluation of r in the parametric region
    if(x < 5.551116e-17)//Linear approximation for r for very small x
      {r = 2*log(x) + 16/3 - log(4) - x;}
    else//Regular evluation for r
      {r = 4.0/(3*a) + 4*a - 4*atanh(a);}
    result[0] = r;
    if(r < -60.0)//10th Order Taylor Series Approximation in parametric region for large negative r
      {fbar = -2*x + 2.5*pow(x,2) - 2.91667*pow(x,3) - 3.28125*pow(x,4) - 3.60938*pow(x,5) - 3.91016*pow(x,6) - 4.18945*pow(x,7) - 4.45129*pow(x,8) - 4.70117*pow(x,9) - 4.93352*pow(x,10);}
    else if(r >= -60.0 && r <= -1.0)//Regular evaluation in parametric region
      {fbar = 4.0/3 - 4.0/(3*pow(a,1.5));}
    else if(r > -1.0 && r < -0.6)//Evaluation in glue region
      {fbar = -13.49319*pow(r,3) - 33.45508*pow(r,2) - 26.57345*r - 6.81477;}
    else// r >= -0.6, in linear region
      {fbar = -r - 0.6;}
    result[1] = fbar;
    return resultPtr;
}

/*Function to calculate a difference quotient for the 4 parameter values
   fbar' approximately = (fbarf - fbar) / (rf - r)
   This function can be used to calculate both the first and second derivatives of
   fbar
*/
double diffQuotient(double r, double rf, double fbar, double fbarf)
{
    return (fbarf - fbar) / (rf - r);
}


//Function required to interface with Python
/*
static PyObject* diff2(PyObject* self, PyObject* args)
{
	return NULL;
}*/

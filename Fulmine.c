/*Eric Gelphman
  University of California San Diego(UCSD)
  Department of Mathematics
  Jacobs School of Engineering Department of Electrical and Computer Engineering(ECE)
  August 29, 2017

  Fulmine, a C subroutine to perform numerical differentiation and evaluation     for fGenerator and KrazyGlue.
*/

//#include <Python.h>//Library needed to interface with Python
#include <stdio.h>
#include <math.h>

/*Function to perform numerical differentiation for first and second order derivatives
  Paramters: r_li = r value at index i-1, r_ui = r value at index i+1, f_li = fbar value
  at index i-1, f_c = fbar value at index i, f_ui = fbar value at index i+1
  Returns a pointer to an array of length 2. fbar' is at index 0, fbar'' is at index 1
*/
double* numericDiff(double r_li, double r_ui, double f_li, double f_c, double f_ui)
{
    double deriv[2];
    double* derivptr;
    derivptr = deriv;
    deriv[0] = (f_ui - f_li) / (2*(r_ui - r_li)); //First derivative fbar'
    deriv[1] = (f_ui - 2.0*f_c + f_li) / ((r_ui - r_li)*(r_ui - r_li));//Second derivative fbar''
    return derivptr;
}

/*Function to evaluate r and fbar in all 3 regions: parametric (r <= -1.0), glue region
(-1.0 < r < -0.6), and linear(r >= -0.6)
*/
double* rfBarEval(double x)
{
    double* rfPtr;
    double rfBar[2];
    rfPtr = rfBar;
    double a,r,fbar;
    a = 1 - x;
    //Evaluation of r in the parametric region
    if(x < 5.551116e-17)//Linear approximation for r for very small x
      {r = 2*log(x) + 16/3 - log(4) - x;}
    else//Regular evluation for r
      {r = 4.0/(3*a) + 4*a - 4*atanh(a);}
    rfBar[0] = r;
    //Evaluation of fbar
    if(r < -60.0)//10th Order Taylor Series Approximation in parametric region for large negative r
      {fbar = -2*x + 2.5*pow(x,2) - 2.91667*pow(x,3) - 3.28125*pow(x,4) - 3.60938*pow(x,5) - 3.91016*pow(x,6) - 4.18945*pow(x,7) - 4.45129*pow(x,8) - 4.70117*pow(x,9) - 4.93352*pow(x,10);}
    else if(r >= -60.0 && r <= -1.0)//Regular evaluation in parametric region
      {fbar = 4.0/3 - 4.0/(3*pow(a,1.5));}
    else if(r > -1.0 && r < -0.6)//Evaluation in glue region
      {fbar = -13.49319*pow(r,3) - 33.45508*pow(r,2) - 26.57345*r - 6.81477;}
    else// r >= -0.6, in linear region
      {fbar = -r - 0.6;}
    rfBar[1] = fbar;
    return rfPtr;
}

//Function required to interface with Python
/*
static PyObject* diff2(PyObject* self, PyObject* args)
{
	return NULL;
}*/

int main()
{
  return 0;
}

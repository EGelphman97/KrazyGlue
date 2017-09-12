/*
Eric Gelphman
University of California San Diego(UCSD)
Department of Mathematics
Irwin and Joan Jacobs School of Engineering Department of Electrical and Computer Engineering(ECE)
September 12, 2017

lightPoly, a C library that does the operations of numerical evaluation and differentiation
on univariate polynomials
*/

#include<stdio.h>
#include<stdlib.h>

//Structure definition to create a polynomial data type
typedef struct polynomial{
  double* coefs;//Array of coefficients in ascending degree order
  int deg;//Degree of polynomial, length of coefs = deg + 1
}poly;

/*
Function that implement's Horner's Rule to evaluate a polynomial at a pointer
Parameters: coefs is an array of the polynomial coefficients, in ascending degree order,
deg is the degree of the polynomial, val is the value the polynomial is to be evaluated at
Return: value of polynomial at point x = val
*/
double eval(double* coefs, int deg, double val)
{
  double y = 0;//Value of polynomial at point val y = f(val)
  int k = deg;
  while(k >= 0)
  {
    y = coefs[k] + val*y;
    k -= 1;
  }
  return y;
}

/*
Function that differentiates a polynomial using the power rule
Parameters: coefs is an array of the polynomial coeffficients in ascending degree order
deg is the degree of the polynomial
Return: array of polynomial coefficients of the derivative in ascending degree order
*/
double* diff(double* coefs, int deg)
{
  double* result = malloc(sizeof(double)*deg);
  int i;
  for(i = 1; i <= deg; i += 1)
  {
    result[i - 1] = i * coefs[i];
  }
  return result;
}

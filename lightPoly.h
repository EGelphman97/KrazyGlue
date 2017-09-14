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

//Function to initialize a polynomial from a pre-initialized array
//Parameters: coefs1 is the pre-initialized array, deg1 is degree of new polynomial
//Returns pointer to newly initialized polynomial
poly* poly_init(double* coefs1, int deg1)
{
    poly* newPoly = malloc(sizeof(poly));
    newPoly->coefs = coefs1;
    newPoly->deg = deg1;
    return newPoly;
}

//Function to initialize a polynomial using a blank array
//Parameter: deg1 is degree of new polynomial
//Returns pointer to newly initialized polynomial, with its array of coefficients set to all 0.0
poly* poly_init_blank(int deg1)
{
    poly* newPoly = malloc(sizeof(poly));
    double* nCoefs = malloc((deg1 + 1)*sizeof(double));//Build array of coefficients
    int i;
    for(i = 0; i <= deg1; i += 1)
    {
      nCoefs[i] = 0.0;
    }
    newPoly->coefs = nCoefs;
    newPoly->deg = deg1;
    return newPoly;
}

/*
Function that implement's Horner's Rule to evaluate a polynomial at a pointer
Parameters: coefs is an array of the polynomial coefficients, in ascending degree order,
deg is the degree of the polynomial, val is the value the polynomial is to be evaluated at
Return: value of polynomial at point x = val
*/
double eval(poly* polynom, double val)
{
  double y = 0;//Value of polynomial at point val y = f(val)
  int k = polynom->deg;
  while(k >= 0)
  {
    y = polynom->coefs[k] + val*y;
    k -= 1;
  }
  return y;
}

//Function to perform polynomial addition/subtraction
//Parameters: first and second are the two polynomials to be added or subtracted, format is first + second or first - second
//Return: Pointer to the resultant polynomial of the addition or subtraction
poly* polyAdd(poly* first, poly* second, char op)
{
  double* nCoefs;//Array of coefficients of new polynomial
  int highestCommonDeg, nDeg, i;//highest common deg = degree for which both polynomials have terms
  int deg1 = first->deg;
  int deg2 = second->deg;
  if(deg1 >= deg2)//degree of first >= degree of second
  {
    nCoefs = malloc(sizeof(double)*deg1);//Build new coefficient array
    highestCommonDeg = deg2;
    nDeg = deg1;
    for(i = highestCommonDeg + 1; i <= deg1; i += 1)//Fill new coefficient array with higher-degree terms not involved in the addition
    {
        nCoefs[i] = first->coefs[i];
    }
  }
  else//degree of first < second
  {
    nCoefs = malloc(sizeof(double)*deg2);//Build new coefficient array
    highestCommonDeg = deg1;
    nDeg = deg2;
    if(op == '+')//If doing addition
    {
      for(i = highestCommonDeg + 1; i <= deg2; i += 1)//Fill new array with higher-degree terms not involved in the addition
      {
        nCoefs[i] = second->coefs[i];
      }
    }
    else//If doing subtraction
    {
      for(i = highestCommonDeg + 1; i <= deg2; i += 1)//Fill new array with opposites of higher-degree terms not involved in the subtraction
      {
        nCoefs[i] = -(second->coefs[i]);
      }
    }
  }
  if(op == '+')//Perform termwise addition
  {
    for(i = 0; i <= highestCommonDeg; i += 1)
    {
      nCoefs[i] = first->coefs[i] + second->coefs[i];
    }
  }
  else//Perform termwise subtraction
  {
    for(i = 0; i <= highestCommonDeg; i += 1)
    {
      nCoefs[i] = first->coefs[i] - second->coefs[i];
    }
  }
  return poly_init(nCoefs, nDeg);//Initialize new polynomial and return a pointer to it
}

/*
Function that differentiates a polynomial using the power rule
Parameter: polynom is the polynomial to be differentiated
Return: array of polynomial coefficients of the derivative in ascending degree order
*/
poly* diff(poly* polynom)
{
  int deg1 = polynom->deg;
  double* nCoefs = malloc(deg1*sizeof(double));//Build array of new polynomial
  int i;
  for(i = 1; i <= deg1; i += 1)//Perform differentiation
  {
    nCoefs[i - 1] = i * polynom->coefs[i];
  }
  return poly_init(nCoefs, deg1 - 1);//Degree is decreased by 1 due to differentiation
}

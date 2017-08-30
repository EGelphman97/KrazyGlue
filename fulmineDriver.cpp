/*Eric Gelphman
  University of California Sabn Diego(UCSD)
  Department of Mathematics
  Irwin and Joan Jacobs School of Engineering, Department of Electrical and Computer Engineering(ECE)
  August 29, 2017

Driver C++ program to test the evaluation of r and fbar as well as the numerical
calculation of derivatives. This driver program prints its necessary data to a file
*/

#include <iostream>
#include <vector>
#include "Fulmine.c"

int main()
{
    double x;
    vector<double*> data;
    x = 1.0e-308; //initial value of x, is the smallest value a double can hold
    while(x < 0.6)
    {
        double vals[4];
        double* rfPtr, derivPtr, newPtr;
        rfPtr = rfBarEval(x);//Evaluate r and fbar
        deriv = numericDiff(1.0, 2.0, 3.0, 4.0, 5.0);//Numerically calculate derivavives
        newPtr = vals;
        vals[0] = *rfPtr;         //r
        vals[1] = *(rfPtr + 1);   //fbar
        vals[2] = *derivPtr;      //fbar'
        vals[3] = *(derivPtr + 1);//fbar''
        data.push_back(newPtr);   //Append block of data to list
        x *= 10;//Step size is scaled exponentially
    }
    return 0;
}

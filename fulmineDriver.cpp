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
    int idx;
    vector<double*> data;
    x = 1.0e-308; //initial value of x, is the smallest value a double can hold
    idx = 0;
    while(x < 0.6)
    {
        double vals[4];
        double* rfPtr, derivPtr, newPtr;
        rfPtr = rfBarEval(x);//Evaluate r and fbar
        newPtr = vals;
        vals[0] = *rfPtr;         //r
        vals[1] = *(rfPtr + 1);   //fbar
        data.push_back(newPtr);   //Append block of data to list
        //Finite Difference differentiation, starts at index 1 in vector data
        if(idx >= 2)
        {
            double* bPtr, fPtr;
            int i;
            i = idx - 1;
            bPtr = data[i - 1];
            fPtr = data[i + 1]
            deriv = numericDiff(*(bPtr), *(fPtr), *(bPtr + 1), *(rfPtr + 1), *(fPtr + 1));
            *(data[i] + 2) = *derivPtr;      //fbar'
            *(data[i] + 3) = *(derivPtr + 1);//fbar''
        }
        x *= 10;//Step size is scaled exponentially
        idx++;
    }
    return 0;
}

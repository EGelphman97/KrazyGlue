/*
Eric Gelphman
University of California San Diego(UCSD)
Department of Mathematics
Irwin and Joan Jacobs School of Engineering Department of Electrical and Computer Engineering(ECE)
February 1, 2018
Version 1.1.0
*/

#include "fulminePlusPlus2.hpp"
using namespace std;

double TIME_STEP = 1.0e-04;

//Function to calculate the first and second derivatives of a function 5 using finite-difference method
vector< array<double,4> > finiteDiff2(vector< array<double,4> > fvals)
{
  int i;
  double h = fvals[1][0] - fvals[0][0];//step size
  //Calculate first derivative
  for(i = 0; i < fvals.size(); i++)
  {
    double fp;
    if(i < 2)//5-point (4th order) left endpoint finite difference
      fp = (-25.0*fvals[i][1] + 48.0*fvals[i+1][1] - 36.0*fvals[i+2][1] + 16.0*fvals[i+3][1] - 3.0*fvals[i+4][1]) / (12.0*h);
    else if(i >= 2 && i < fvals.size() - 2)//5-point (4th order) midpoint finite difference
      fp = (fvals[i-2][1] - 8.0*fvals[i-1][1] + 8.0*fvals[i+1][1] - fvals[i+2][1]) / (12.0*h);
    else//5-point (4th order) right endpoint finite difference
      fp = (-25.0*fvals[i][1] + 48.0*fvals[i-1][1] - 36.0*fvals[i-2][1] + 16.0*fvals[i-3][1] - 3.0*fvals[i-4][1]) / (12.0*h);
    fvals[i][2] = fp;
  }
  for(i = 0; i < fvals.size(); i++)
  {
    double fpp;
    if(i < 2)//5-point (4th order) left endpoint finite difference
      fpp = (-25.0*fvals[i][2] + 48.0*fvals[i+1][2] - 36.0*fvals[i+2][2] + 16.0*fvals[i+3][2] - 3.0*fvals[i+4][2]) / (12.0*h);
    else if(i >= 2 && i < fvals.size() - 2)//5-point (4th order) midpoint finite difference
      fpp = (fvals[i-2][2] - 8.0*fvals[i-1][2] + 8.0*fvals[i+1][2] - fvals[i+2][2]) / (12.0*h);
    else//5-point (4th order) right endpoint finite difference
      fpp = (-25.0*fvals[i][2] + 48.0*fvals[i-1][2] - 36.0*fvals[i-2][2] + 16.0*fvals[i-3][2] - 3.0*fvals[i-4][2]) / (12.0*h);
    fvals[i][3] = fpp;
  }
  return fvals;
}

/*Function to perform the evolution for each new evolved f
Paramters: fvals = previous f (as set of 4-tuples (z, w = f(z), f'(z), f''(z))) that was evolved
           Hvals = values of (z, H) calculated using previous f Values
Return: vector of arrays of size 4 representing the 4-tuples (z, f(z) = w, f'(z), f''(z)) of the newly evolved f
*/
vector< array<double,4> > evolve_step(vector< array<double,4> > fvals, vector< array<double,2> > Hvals)
{
  vector< array<double,4> > f_next;//evolved f
  int i;
  for(i = 0; i < fvals.size(); i++)
  {
    double w_next;//What is being calculated, the next evolved value of f(z)
    array<double,4> f_cur_element = fvals[i];
    double z_cur = f_cur_element[0];//current z value
    double f_cur = f_cur_element[1];//current value of f
    double fp_cur = f_cur_element[2];//f'(z) at z = z_cur
    double r = 2*LambertW(((z_cur*z_cur) - (f_cur*f_cur)) / exp(1.0)) + 2.0;
    double alpha = (32.0 / r)*exp(-r/2.0);
    double df_ds = (alpha/sqrt(1.0 - (fp_cur*fp_cur))) - sqrt(((alpha*alpha)/(1 - (fp_cur*fp_cur))) - 1.0);
    w_next = f_cur + (TIME_STEP*Hvals[i][1]*df_ds);//step for Euler's method
    array<double,4> f_next_element = {z_cur, w_next, 0.0, 0.0};
    f_next.push_back(f_next_element);
  }
  f_next = finiteDiff2(f_next);
  return f_next;
}

/*Function to perform the evolution
  Parameters: num_evol = number of evolutions to be performed
  Return: final evolved f after num_evol evolutions
*/
vector< array<double,4> > evolve(int num_evol)
{
  int i = 0;
  vector< array<double,4> > evolved_f;//values of (z, f(z) = w, f'(z), f''(z))
  vector< array<double,2> > Hvals;//values of (z, H)
  for(i = 0; i < num_evol; i++)
  {
    if(i == 0)//First step is to do fGenerator to establish initial approximation
    {
      evolved_f = fGeneratorPlusPlus(0.0005, 5000);
      Hvals = calculateH(evolved_f);
    }
    else//General case
    {
      Hvals = calculateH(evolved_f);
      evolved_f = evolve_step(evolved_f, Hvals);
    }
  }
  return evolved_f;
}

int main()
{
  return 0;
}

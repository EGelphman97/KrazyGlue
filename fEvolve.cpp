/*
Eric Gelphman
University of California San Diego(UCSD)
Department of Mathematics
Irwin and Joan Jacobs School of Engineering Department of Electrical and Computer Engineering(ECE)
January 29, 2018
Version 1.0.0
*/

#include "fulminePlusPlus.cpp"
using namespace std;

//Function to computer Euler's method
vector< array<double,2> > euler(double z_start, double z_end, double w_init, double h)
{
  double z, z_prev, w_i, w_i_1;
  int firstflag = 1;
  vector< array<double,2> > sol;//Vector that holds the ordered pairs (z,w) that represent the numerical solution
  z = z_start;
  while(z <= z_end)
  {
    if(firstflag)
    {
      w_i = w_init;
      firstflag = 0;
    }
    else
      w_i = w_i_1 + h*(w_i_1 - (z_prev*z_prev) + 1.0);
    array<double,2> op = {z, w_i};
    sol.push_back(op);
    w_i_1 = w_i;//Set w_(i-1)
    z_prev = z;
    z += h;//Increment z
  }
  return sol;
}

vector< array<double,2> > f_evolve(vector< array<double,4> > fGen, vector< array<double,2> > hVals)
{
  vector< array<double,2> > sol;
  //Check that the size of the parameter vectors are equal
  if(fGen.size() != hVals.size())
  {
    cout << "Fatal error! Size mismatch on parameter vectors from fGenerator and calculation of H!\n";
    exit(1);
  }
  double z, H, z_prev, w_sol, w_sol_prev;
  double step_size = fGen[1][0] - fGen[0][0];
  double w_sol_init = 0.0;//Initial value for solution
  int i;
  for(i = 0; i < fGen.size(); i++)//Iterate through vector
  {
    if(i == 0)//Special case for first element
    {
      z = fGen[0][0];
      w_sol = w_sol_init;
    }
    else//General case-apply Euler's method to solve the IVP
    {
      array<double,4> fGen_element = fGen[i];
      array<double,2> hVals_element = hVals[i];
      z = fGen_element[0];
      H = hVals_element[0];//value of H at current z-value
      w_sol = w_sol_prev + step_size*H;
    }
    array<double,2> sol_element = {z, w_sol};
    sol.push_back(sol_element);
    w_sol_prev = w_sol;
  }
  return sol;
}

int main()
{
  return 0;
}

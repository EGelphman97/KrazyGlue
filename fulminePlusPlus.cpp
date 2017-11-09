/*Eric Gelphman
  University of California, San Diego(UCSD)
  Department of Mathematics
  Irwin and Joan Jacobs School of Engineering Department of Electrical and Computer Engineering(ECE)

  November 9, 2017
  fulminePlusPlus, a C++ extension that performs all mathematical and file I/0 operations for KG2

  Version 1.0.0
*/

#include <iostream>
#include <vector>
#include <cmath>
using namespace std;

double ALPHA = 0.8;//rho parameters
double BETA = 1.2;
double GAMMA = 2.0;

//Function to generate the 4-tuple (z, f(z) = w, f'(z), f''(z)) in the Schwarzschild region
vector<double> genfS(double z)
{
  vector<double> result;
  result.push_back(z);
  result.push_back(0.0);//f(z)
  result.push_back(0.0);//f'(z)
  result.push_back(0.0);//f''(z)
  return result;
}

//Function to generate the 4-tuple (z, f(z) = w, f'(z), f''(z)) in the Glue Region
vector<double> genfG(double z)
{
  vector<double> result;
  return result;
}

//Function to generate the 4-tuple (z, f(z) = w, f'(z), f''(z)) in the Friedman Region
vector<double> genfF(double z)
{
  vector<double> result;
  return result;
}

/*Function to generate a vector of 4-tuples (z, f(z) = w, f'(z), f''(z)) for allocate
3 regions. The 4-tuples are stored horizontally in the vector(some modulus operations are going to
be needed in Python)*/
vector<double> fGeneratorPlusPlus()
{
  double z = 0.0;
  vector<double> result;
  while(z <= GAMMA)
  {
    vector<double> tuple4;
    if(z <= ALPHA)//Schwarzschild region
      tuple4 = genfS(z);
    else if(z > ALPHA && z < BETA)//Glue region
      tuple4 = genfG(z);
    else if(z >= BETA)//Friedman region
      tuple4 = genfF(z);
    result.push_back(tuple4[0]);//Store 4-tuples horizontally
    result.push_back(tuple4[1]);
    result.push_back(tuple4[2]);
    result.push_back(tuple4[3]);
  }
  //Output to file
  return result;//Return list of 4-tuples(to Python)
}

//Function to transform an input 4-tuple (x1, f(x1), f'(x1), f''(x2)) to a 4-tuple (x2, f(x2), f'(x2), f''(x2)) in different Coordinates
//The variables x1 and x2 represent rho(Tolman-Bondi coordinates) and z(Kruskal coordinates)
vector<double> transformCoordinates(vector<double> tuple4, char domainCoords, char codomainCoords)
{

  if(domainCoords == 'T' && codomainCords == 'K')
  {
    //Transform the input 4-tuple from Tolman-Bondi coordinates to Kruskal coordinates
  }
  else if(domainCoords == 'K' && codomainCords == 'T')
  {
    //Transfrom the input 4-tuple from Kruskal Coordinates to Tolman-Bondi coordinates
  }
  else//Invalid coordinate change
  {
    cout << "Error! Invalid Coordinates, Only transformations Tolman-Bondi <==> Kruskal are allowed.\n";
    return NULL;
  }
  return tuple4;
}

//Function to calculate the mean curvature H(z,w) in the Schwarzschild region
double calcHS(double z, double f, double fp, double fpp)
{
  return 0.0;
}

//Function to calculate the mean curvature H(z,w) in the Glue region
double calcHG(double z, double f, double fp, double fpp)
{
  return 0.0;
}

//Function to calculate the mean curvature H(z,w) in the Friedman region
double calcHF(double z, double f, double fp, double fpp)
{
  return 0.0;
}

//Function to calculate the mean curvature H for all three regions
//Returns a vector that stores the values of H as well as their corresponding z-values
vector<double> calculateH(vector<double> fData)
{
  vector<double> h_zValues;
  return h_zValues;
}

int main()
{
  return 0;
}

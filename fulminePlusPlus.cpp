/*Eric Gelphman
  University of California, San Diego(UCSD)
  Department of Mathematics
  Irwin and Joan Jacobs School of Engineering Department of Electrical and Computer Engineering(ECE)
  January 15, 2018
  fulminePlusPlus, a C++ extension that performs all mathematical and file I/0 operations for KG2
  Version 2.0.1
*/

#include <stdio.h>
#include <iostream>
#include <vector>//std::vector
#include <array>//std::array
#include <cmath>
using namespace std;

double BETA = 0.75;//End of gluing region for g(rho)
double METRIC_GLUE_END_TB = 1.25;//Right endpoint of Glue region for the metric
double GAMMA = 2.0;
double A = 0.40;//z parameters
double B = 1.691486380082176;//When rho = BETA, z = B
double METRIC_GLUE_START_K = B;//Start of gluing of metric in Kruskal coordinates
double EPSILON = 1.0e-15;//Slightly better than machine precision
double TAUFRIEDMAN = -1.7;//Tau in region where tau is constant

//Function to compute the lambert W-function
double LambertW(const double z);
const int dbgW=0;
double LambertW(const double z) {
  int i;
  const double eps=4.0e-16, em1=0.3678794411714423215955237701614608;
  double p,e,t,w;
  if (dbgW) fprintf(stderr,"LambertW: z=%g\n",z);
  if (z<-em1 || isinf(z) || isnan(z)) {
    fprintf(stderr,"LambertW: bad argument %g, exiting.\n",z); exit(1);
  }
  if (0.0==z) return 0.0;
  if (z<-em1+1e-4) { // series near -em1 in sqrt(q)
    double q=z+em1,r=sqrt(q),q2=q*q,q3=q2*q;
    return
     -1.0
     +2.331643981597124203363536062168*r
     -1.812187885639363490240191647568*q
     +1.936631114492359755363277457668*r*q
     -2.353551201881614516821543561516*q2
     +3.066858901050631912893148922704*r*q2
     -4.175335600258177138854984177460*q3
     +5.858023729874774148815053846119*r*q3
     -8.401032217523977370984161688514*q3*q;  // error approx 1e-16
  }
  /* initial approx for iteration... */
  if (z<1.0) { /* series near 0 */
    p=sqrt(2.0*(2.7182818284590452353602874713526625*z+1.0));
    w=-1.0+p*(1.0+p*(-0.333333333333333333333+p*0.152777777777777777777777));
  } else
    w=log(z); /* asymptotic */
  if (z>3.0) w-=log(w); /* useful? */
  for (i=0; i<10; i++) { /* Halley iteration */
    e=exp(w);
    t=w*e-z;
    p=w+1.0;
    t/=e*p-0.5*(p+1.0)*t/p;
    w-=t;
    if (fabs(t)<eps*(1.0+fabs(w))) return w; /* rel-abs error */
  }
  /* should never get here */
  fprintf(stderr,"LambertW: No convergence at z=%g, exiting.\n",z);
  exit(1);
}

/*Helper function to tranfrom an input ordered pair (z, f(z) = w) in Kruskal coordinates to
an output 4-tuple (rho, g(rho) = tau) in Tolman-Bondi Coordinates*/
array<double,2> KTBTransform2(array<double,2> vals)
{
  double z = vals[0];
  double w = vals[1];//f(z)
  double U = w - z;//U = w - z
  double V = w + z;//V = w + z
  double t = 2.0*log(-V / U);
  double r = 2.0*LambertW((-U*V) / exp(1)) + 2.0;
  double tau = 2.0*sqrt(2.0)*sqrt(r) + t - 2.0*log(fabs((sqrt(2.0) + sqrt(r))/(sqrt(2.0) - sqrt(r))));
  double rho = tau + sqrt((2.0/9.0)*pow(r,3));
  array<double,2> result = {rho, tau};
  return result;
}

/*Function to compute the transformation T: (rho, g(rho) = tau) ---> (z, f(z) = w),
where each ordered pair is represented by an array of size 2
*/
array<double,2> discTBKTransform(array<double,2> rhotau)
{
  double z, w;//What is going to be returned
  double rho = rhotau[0];
  double tau = rhotau[1];
  double r = pow(4.5, 1.0/3.0)*pow(rho - tau, 2.0/3.0);
  double sqr = sqrt(r);
  double a = r + 2.0*log((r/2.0) - 1);
  double e_t_4 = exp(-(sqr/sqrt(2.0)) + tau/4.0)*sqrt(fabs((sqrt(2.0) + sqr)/(sqrt(2.0) - sqr)));
  double U = (-exp(a/4.0)) / e_t_4;
  double V = exp(a/4.0) * e_t_4;
  z = (V - U) / 2.0;
  w = (V + U) / 2.0;
  array<double,2> result = {z, w};//Build output ordered pair
  return result;
}

//Recursive helper function for binarySearch()
array<double,2> binSearch(int low, int high, vector< array<double,2> > vals, double key)
{
  if(low > high)//Exact Match not found, return indices for linear interpolation
  {
    int idx1, idx2;
    if(low == high)
    {
      idx1 = high;
      idx2 = high + 1;
    }
    else//low > high
    {
      idx1 = high;
      idx2 = low;
    }
    array<double,2> result = {idx1, idx2};
    return result;
  }
  int mid = (high + low) / 2;
  if(key < vals[mid][0])//key < value at mid
    return binSearch(low, mid - 1, vals, key);
  else if(fabs(key - vals[mid][0]) < EPSILON)//Exact Match
  {
    array<double,2> result = {mid, mid};//Average of two equal numbers is that number
    return result;
  }
  else//key > value at mid
    return binSearch(mid + 1, high, vals, key);
}

/*Function to do the binary search, where key is some z value
If key is found, return array of size two, with both indices having the index where key is found
If key is not found, return array of size 2 <id1, id2> , from which to construct a linear approximation*/
array<double,2> binarySearch(vector< array<double,2> > vals, double key)
{
  array<double,2> result;
  result = binSearch(0, vals.size() - 1, vals, key);
  return result;
}
/*Function to search through the vector vals to determine a value for f(z) = w at
the z-value z = target via a linear interpolation of the values for w at the two z-Values
closest to z0
*/
double interpolate(double key, vector< array<double,2> > vals)
{
  double val;
  if(key > vals[vals.size() - 1][0])//Special case
    val = vals[vals.size() - 1][1];
  else//Linear interpolation
  {
    array<double,2> indices = binarySearch(vals, key);
    double x1 = vals[indices[0]][0];
    double y1 = vals[indices[0]][1];
    double m = (vals[indices[1]][1] - y1) / (vals[indices[1]][0] - x1);//Calculate slope
    double B1 = -m*x1 + y1;//Calculate y-intercept
    val = m*key + B1;//y = mx + b
  }
  return val;
}

/*Function to generate a vector of ordered pairs(implemented as arrays of size 2) (z, w)
with an even z-step given by the value of the variable z_step. Variable nPoints
is needed for the initial generation of rho-values*/
vector< array<double,2> > getConstZStep(double z_step, int nPoints)
{
  vector< array<double,2> > result;
  //Generate values of rho and tau for ALPHA < rho <= GAMMA
  double rho_step = (GAMMA - BETA) / (double)nPoints;
  double rho = BETA;
  double tau;
  vector< array<double,2> > zw;
  int i;
  for(i = 0; i < nPoints; i++)//Iterate based on list index i, increment rho seperately
  {
    tau = TAUFRIEDMAN;
    array<double,2> rhotau = {rho, tau};// (rho, tau) to be transformed
    array<double,2> element = discTBKTransform(rhotau);//What is going to be returned
    zw.push_back(element);//Append array to vector
    rho += rho_step;
  }
  double z = zw[0][0];
  array<double,2> rhotauEnd = {GAMMA, TAUFRIEDMAN};
  array<double,2> zwEnd = discTBKTransform(rhotauEnd);
  double z_end = zwEnd[0];//z value to signal the end of the iteration
  int firstFlag = 1;
  while(z < z_end)
  {
    array<double,2> element;
    if(firstFlag == 1)//First element is special case
    {
      element[0] = z;
      element[1] = zw[0][1];//Corresponding w value for z0 = zw[0][0]
      result.push_back(element);//z, element at index 0 in zw
      firstFlag = 0;//Case handled
    }
    else//General case
    {
      double w = interpolate(z, zw);
      element[0] = z;
      element[1] = w;
      result.push_back(element);
    }
    z += z_step;
  }
  return result;
}

/*Function to compute the transfromation T: (rho, g(rho) = tau, g'(rho), g''(rho)) |--> (z, f(z) = w, f'(z), f''(z))
Returns a vector of arrays of size 4, with the array of size 4 representing the 4-tuple (z, f(z) = w, f'(z), f''(z))
Parameters: step_size is the constant z-step, nPoints is a paramter necessary to generate values of (rho, tau)*/
vector< array<double,4> > TBKTransform(double z_step, int nPoints)
{
  vector< array<double,2> > constZStep = getConstZStep(z_step, nPoints);//Transform (rho, g(rho)) |--> (z, f(z))
  vector< array<double,4> > result;
  //Transform (g'(rho), g''(rho)) |--> (f'(z), f''(z))
  int i;
  double h = z_step;//step size
  for(i = 0; i < constZStep.size(); i++)
  {
    array<double,4> element;
    double z = constZStep[i][0];
    double w = constZStep[i][1];
    double fp, fpp;//f'(z), f''(z)
    if(i < 2)//First 2 points
    {
      fp = (-3.0*w + 4.0*constZStep[i + 1][1] - constZStep[i + 2][1]) / (2.0*h);//3 - point endpoint
      double fpplus1 = (constZStep[i + 2][1] - w) / (2.0*h);
      double fpplus2 = (constZStep[i + 3][1] - constZStep[i + 1][1]) / (2.0*h);
      fpp = (-3.0*fp + 4.0*fpplus1 - fpplus2) / (2.0*h);//Compute f'' in similar manner to f'
    }
    else if(i >= constZStep.size() - 2)//Last 2 points
    {
      fp = (-3.0*w + 4.0*constZStep[i - 1][1] - constZStep[i - 2][1]) / (2.0*h);//3 - point endpoint
      double fpminus1 = (w - constZStep[i - 2][1]) / (2.0*h);
      double fpminus2 = (constZStep[i - 3][1] - constZStep[i - 1][1]) / (2.0*h);
      fpp = (-3.0*fp + 4.0*fpminus1 - fpminus2) / (2.0*h);//Compute f'' in similar manner to f'
    }
    else//Use midpoint formulas
    {
      fp = (constZStep[i - 2][1] - 8.0*constZStep[i - 1][1] + 8.0*constZStep[i + 1][1] - constZStep[i + 2][1]) / (12.0*h);//5-point midpoint
      fpp = (constZStep[i - 1][1] - 2.0*w + constZStep[i + 1][1]) / (h*h);//3-point midpoint
    }
    element[0] = z;
    element[1] = w;
    element[2] = fp;
    element[3] = fpp;
    result.push_back(element);
  }
  return result;
}

/*Function to generate a vector of double arrays of size 4 representing values of z, w, and the first and second derivatives.
Each array represents the 4-tuple (z, f(z) = w, f'(z), f''(z)).
Parameters: step_size is the constant z-step, nPoints is a parameter needed to first generate ordered pairs (rho, tau) to start the
interpolation*/
vector< array<double,4> > fGeneratorPlusPlus(double step_size, int nPoints)
{
  vector< array<double,4> > result;
  double z;
  //Generate values for 0.0 <= z <= B
  for(z = 0.0; z <= B; z += step_size)
  {
      if(z <= A)
      {
        array<double, 4> element = {z, 0.0, 0.0, 0.0};
        result.push_back(element);
      }
      else//Glue region for f(z)
      {
        double f, fp, fpp;
        f = 0.12636326487700675 - 1.0428023341486672*z + 3.102608302347679*pow(z,2) - 3.9087348740648844*pow(z,3) + 1.8645894143067527*pow(z,4) - 0.31583170500671803*pow(z,5);
        fp = -1.0428023341486672 + 6.205216604695358*z - 11.726204622194654*pow(z,2) + 7.458357657227011*pow(z,3) - 1.5791585250335902*pow(z,4);
        fpp = 6.205216604695358 - 23.452409244389308*z + 22.375072971681032*pow(z,2) - 6.316634100134361*pow(z,3);
        array<double, 4> elem = {z, f, fp, fpp};
        result.push_back(elem);
      }
  }
  //Generate values for B < z <= C, where C is the corresponding value of z for rho = GAMMA
  vector< array<double,4> > fGlueS = TBKTransform(step_size, nPoints);
  int i;
  for(i = 0; i < fGlueS.size(); i++)//Append contents of fGlueS to result
  {
    result.push_back(fGlueS[i]);
  }
  return result;//Return vector
}

//Function to calculate the mean curvature H(z,w) in the Schwarzschild region
double calcHS(array<double,4> fVals)
{
  double z = fVals[0];
  double f = fVals[1];
  double fp = fVals[2];
  double fpp = fVals[3];
  double r = 2.0*LambertW(((z*z) - (f*f)) / exp(1)) + 2.0;
  double H_Num = exp(r/4.0)*(((r-6.0)*(1.0 - (fp*fp))*(exp(-r/2.0)/(r-1))*(fp*z - f)) + r*fpp);
  double H_Denom = 4*sqrt(2.0)*pow(1.0 - (fp*fp), 3.0/2.0)*sqrt(r);
  double H = H_Num / H_Denom;
  return H;
}

/*Function to calculate the mean curvature H in the Gluing and Friedman region from
  4-tuple (z, f(z) = w, f'(z), f''(z)) represented by kcoords*/
double calcHGF(array<double,4> kcoords)
{
  array<double,2> zw = {kcoords[0], kcoords[1]};
  array<double,2> rhotau = KTBTransform2(zw);
  double r = rhotau[0];
  double t = TAUFRIEDMAN;
  double fp = 0.0;
  double fpp = 0.0;
  double X, Y, dX_dr, dX_dt, dY_dr, dY_dt;
  static double H;
  if(r >= METRIC_GLUE_END_TB)//Friedman Region M(r) = r^3, t0(r) = 1, applies for all r >= 1.2
  {
    X = (pow(3,0.6666666666666666)*pow(r,2)*(1 - t))/(pow(2,0.3333333333333333)*pow(pow(r,6)*(1.0 - t),0.3333333333333333));
    dX_dr = 0.0;
    dX_dt = -((pow(2,0.6666666666666666)*pow(r,2))/(pow(3,0.3333333333333333)*pow(-(pow(r,6)*(-1.0 + t)),0.3333333333333333)));
    Y = pow(4.50, 1.0/3.0)*r*pow(1 - t, 2.0/3);
    dY_dr = pow(4.50, 1.0/3.0)*pow(1 - t, 2.0/3);
    dY_dt = (-1.100642416*r) / pow(1 - t, 1.0/3);
  }
  else//Glue Region applies for 0.75 < r < 1.25
  {
    X = (-35.61435406495837 + 11913.766254916665*pow(r,5) - 7324.775280463868*pow(r,6) + 2450.03001867957*pow(r,7) - 343.4004338849967*pow(r,8) - 3.111571796362811e-12*pow(r,9) + pow(r,3)*(6637.002751810338 - 107.3126355890683*t) + r*(425.85671576491285 - 47.981130335493084*t) +
        5.804168992195918*t + pow(r,4)*(-11456.655701397658 + 33.019272488944544*t) + pow(r,2)*(-2264.722249988459 + 115.56745371130128*t))/pow(pow(1.3427734374994298 - 10.54687499999697*r + 43.593749999993676*pow(r,2) - 69.99999999999355*pow(r,3) + 48.74999999999679*pow(r,4) -
        11.999999999999375*pow(r,5),2)*(2.9531249999999623 - 12.49999999999981*r + 22.49999999999962*pow(r,2) - 15.999999999999632*pow(r,3) + 3.9999999999998277*pow(r,4) + 3.141167326248251e-14*pow(r,5) - t),0.3333333333333333);
    dX_dr = (3020.927028246527 - 3.1112245727108693e8*pow(r,19) + 3.6259915968088605e7*pow(r,20) - 1.9779864991767928e6*pow(r,21) - 3.727911937785611e-8*pow(r,22) - 1.6889416116496277e-22*pow(r,23) + pow(r,14)*(6.015040393797038e10 - 2.4672565858776075e8*t) +
            pow(r,16)*(1.703713323419403e10 - 9.36785817470958e6*t) + pow(r,18)*(1.6605211766910415e9 + 7.169061418819169e-9*t) - 1444.9044797687486*t + 95.13774352388168*pow(t,2) + pow(r,17)*(-6.179564652060596e9 + 654257.0728047623*t) + pow(r,15)*(-3.6093215640140816e10 + 6.1594547132135525e7*t) +
            pow(r,7)*(-5.574060233535852e9 + 6.412153131929202e8*t - 2.3025902607623823e6*pow(t,2)) + pow(r,9)*(-3.092273938222233e10 + 1.8756433194791079e9*t - 2.269553898220539e6*pow(t,2)) + pow(r,5)*(-4.205554691678066e8 + 7.814083855654563e7*t - 584384.8951669121*pow(t,2)) +
            pow(r,11)*(-7.508122754287083e10 + 1.9505625291768303e9*t - 469212.11688596057*pow(t,2)) + pow(r,3)*(-1.1807106966068413e7 + 3.215770737640476e6*t - 49129.63594292593*pow(t,2)) + pow(r,13)*(-8.00852747171778e10 + 6.734023940549386e8*t - 9509.550476815035*pow(t,2)) +
            r*(-86657.1392636461 + 33836.2562642282*t - 1317.9887118348179*pow(t,2)) + pow(r,2)*(1.250059744850459e6 - 406866.0548882588*t + 9777.575055275804*pow(t,2)) + pow(r,12)*(8.60501790946082e10 - 1.3268132050906591e9*t + 100444.62691135757*pow(t,2)) + pow(r,4)*(8.075348292716156e7 - 1.830107320499192e7*t +
            189570.60206377413*pow(t,2)) + pow(r,10)*(5.336749609824144e10 - 2.1812858579780464e9*t + 1.2808284142947402e6*pow(t,2)) + pow(r,6)*(1.7154556879061031e9 - 2.5501246906960225e8*t + 1.3658524645672236e6*pow(t,2)) + pow(r,8)*(1.4582176144399971e10 - 1.2465484036576378e9*t +
            2.73915228255279e6*pow(t,2)))/(3.*pow(pow(1.3427734374994298 - 10.54687499999697*r + 43.593749999993676*pow(r,2) - 69.99999999999355*pow(r,3) + 48.74999999999679*pow(r,4) - 11.999999999999375*pow(r,5),2)*(2.9531249999999623 - 12.49999999999981*r + 22.49999999999962*pow(r,2) -
            15.999999999999632*pow(r,3) + 3.9999999999998277*pow(r,4) + 3.141167326248251e-14*pow(r,5) - t),1.3333333333333333));
    dX_dt = (47.999999999995*pow(0.11189778645829164 - 0.8789062499997933*r + 3.6328124999996625*pow(r,2) - 5.8333333333330994*pow(r,3) + 4.062499999999944*pow(r,4) - 1.*pow(r,5),2)*(15.80695560027668 - 2691.0707078489486*pow(r,5) + 1441.8415653505617*pow(r,6) - 422.6466878584847*pow(r,7) +
            52.830835982310646*pow(r,8) - 4.0389678347315804e-28*pow(r,9) + pow(r,2)*(950.2069547892688 - 231.1349074226026*t) + pow(r,4)*(3033.645659921751 - 66.03854497788907*t) - 11.608337984391834*t + r*(-216.88244800843432 + 95.96226067098618*t) + pow(r,3)*(-2164.826052556411 +
            214.6252711781366*t)))/pow(pow(1.3427734374994298 - 10.54687499999697*r + 43.593749999993676*pow(r,2) - 69.99999999999355*pow(r,3) + 48.74999999999679*pow(r,4) - 11.999999999999375*pow(r,5),2)*(2.9531249999999623 - 12.49999999999981*r + 22.49999999999962*pow(r,2) -
            15.999999999999632*pow(r,3) + 3.9999999999998277*pow(r,4) + 3.141167326248251e-14*pow(r,5) - t),1.3333333333333333);
    Y = 1.6509636244473134*pow((1.3427734374994298 - 10.54687499999697*r + 43.593749999993676*pow(r,2) - 69.99999999999355*pow(r,3) + 48.74999999999679*pow(r,4) - 11.999999999999375*pow(r,5))*pow(2.9531249999999623 - 12.49999999999981*r + 22.49999999999962*pow(r,2) - 15.999999999999632*pow(r,3) +
        3.9999999999998277*pow(r,4) + 3.141167326248251e-14*pow(r,5) - 1.*t,2),0.3333333333333333);
    dY_dr = (-3.111571796362811e-12*(2.9531249999999623 - 12.49999999999981*r + 22.49999999999962*pow(r,2) - 15.999999999999632*pow(r,3) + 3.9999999999998277*pow(r,4) + 3.141167326248251e-14*pow(r,5) - 1.*t)*(1.1445776088660021e13 - 3.828857900320007e15*pow(r,5) + 2.3540434737922385e15*pow(r,6) -
            7.873930537432776e14*pow(r,7) + 1.1036236871873098e14*pow(r,8) + 1.*pow(r,9) + pow(r,2)*(7.278386610380471e14 - 3.714118178034355e13*t) + pow(r,4)*(3.6819512616709055e15 - 1.0611766222955725e13*t) - 1.865349531378497e12*t + r*(-1.368622495751847e14 + 1.5420222792731104e13*t) +
            pow(r,3)*(-2.133006463025692e15 + 3.448824022460563e13*t)))/pow((1.3427734374994298 - 10.54687499999697*r + 43.593749999993676*pow(r,2) - 69.99999999999355*pow(r,3) + 48.74999999999679*pow(r,4) - 11.999999999999375*pow(r,5))*pow(2.9531249999999623 - 12.49999999999981*r + 22.49999999999962*pow(r,2) -
            15.999999999999632*pow(r,3) + 3.9999999999998277*pow(r,4) + 3.141167326248251e-14*pow(r,5) - 1.*t,2),0.6666666666666666);
    dY_dt = (13.20770899557782*(-0.11189778645829164 + 0.8789062499997933*r - 3.6328124999996625*pow(r,2) + 5.8333333333330994*pow(r,3) - 4.062499999999944*pow(r,4) + 1.*pow(r,5))*(2.9531249999999623 - 12.49999999999981*r + 22.49999999999962*pow(r,2) - 15.999999999999632*pow(r,3) +
            3.9999999999998277*pow(r,4) + 3.141167326248251e-14*pow(r,5) - 1.*t))/pow((1.3427734374994298 - 10.54687499999697*r + 43.593749999993676*pow(r,2) - 69.99999999999355*pow(r,3) + 48.74999999999679*pow(r,4) - 11.999999999999375*pow(r,5))*pow(2.9531249999999623 -
            12.49999999999981*r + 22.49999999999962*pow(r,2) - 15.999999999999632*pow(r,3) + 3.9999999999998277*pow(r,4) + 3.141167326248251e-14*pow(r,5) - 1.*t,2),0.6666666666666666);
  }
  double H_Denom = X * Y * pow(pow(X,2) - pow(fp,2), 1.5);
  if(isnan(H_Denom))//If there is a negative number under the radical, throw an exception
  {
      printf("Error! Cannot Have a Negative Number Raised to the Power 3/2!\n");
      printf("r= %lf f'(r)= %lf X= %lf\n", r, fp, X);
      exit(1);
  }
  else if(H_Denom == 0.0)//If the denominator = 0, throw an exception
  {
      printf("Error! Cannot Divide by 0!\n");
      printf("r= %lf f'= %lf X= %lf\n", r, fp, X);
      exit(1);
  }
  double H_Num = (2.0*pow(fp,3)*dY_dr) - (pow(X,3)*Y*dX_dt) + ((X*Y*fp)*(dX_dr + 2*fp*dX_dt)) - (2.0*pow(X,4)*dY_dt) -(pow(X,2)*((Y*fpp)+(2.0*fp)*(dY_dr - fp*dY_dt)));
  H = H_Num / H_Denom;
  return H;
}

/*Function to calculate the mean curvature H for all three regions
Returns a vector that stores the values of H as well as their corresponding z-values
The ordered pair (z, H) is represented by an array of size 2*/
vector< array<double,2> > calculateH(vector< array<double,4> > fVals)
{
  vector< array<double,2> > hValues;
  int i;
  double H;
  for(i = 0; i < fVals.size(); i++)
  {
    array<double,4> element = fVals[i];
    double z = element[0];
    array<double,2> zw = {z, element[1]};
    array<double,2> rt = KTBTransform2(zw);
    if(z <= METRIC_GLUE_START_K)//Schwarzschild region, rho <= 0.75
    {
      H = calcHS(element);
    }
    else//Glue and Friedman Regions
    {
      H = calcHGF(element);
    }
    //printf("rho: %lf z: %lf H: %lf\n", rt[0], z, H);
    array<double,2> nElem = {element[0], H};//Build ordered pair
    hValues.push_back(nElem);//Append ordered pair to vector
  }
  return hValues;
}

/*Function to output the numerical data from fGeneratorPlusPlus to a file*/
void outputToFileF(vector< array<double,4> > fGen)
{
  FILE* file1;
  file1 = fopen("fGenerator5.txt", "w");
  int i;
  for(i = 0; i < fGen.size(); i++)
  {
      array<double,4> element = fGen[i];
      /*Print formatted string in this format: z f(z) = w f'(z) f''(z) where each
      value is separated by a space*/
      fprintf(file1, "%.15lf %.15lf %.15lf %.15lf\n", element[0], element[1], element[2], element[3]);
  }
  fclose(file1);
}

/*Function to output the numerical data from calculateH to a file*/
void outputToFileH(vector< array<double,2> > h_values)
{
  FILE* file1;
  file1 = fopen("outputH.txt", "w");
  int i;
  for(i = 0; i < h_values.size(); i++)
  {
      array<double,2> element = h_values[i];
      /*Print formatted string in this format: z H where z and H are separated by a space*/
      fprintf(file1, "%.15lf %.15lf\n", element[0], element[1]);
  }
  fclose(file1);
}

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

int main()
{
  //vector< array<double,4> > fGen = fGeneratorPlusPlus(0.0005, 5000);//fGenerator's Newest Form
  //outputToFileF(fGen);//Output 4-tuples (z, f(z) = w, f'(z), f''(z)) to file
  //vector< array<double,2> > h_vals = calculateH(fGen);
  //outputToFileH(h_vals);//Output ordered pairs (z, H) to file
  vector< array<double,2> > sol1 = euler(0.0, 2.0, 0.5, 0.5);
  int i;
  printf("Numerical solution of dw/dz = w - z^2 + 1\n");
  for(i = 0; i < sol1.size(); i++)
  {
    printf("z: %lf w: %lf\n", sol1[i][0], sol1[i][1]);
  }
  return 0;
}

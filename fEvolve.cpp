/*
Eric Gelphman
University of California San Diego(UCSD)
Department of Mathematics
Irwin and Joan Jacobs School of Engineering Department of Electrical and Computer Engineering(ECE)
February 3, 2018
Version 1.2.0
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
    double fp, fpp;
    if(i < 2)//5-point (4th order) left endpoint finite difference
      fp = (-25.0*fvals[i][1] + 48.0*fvals[i+1][1] - 36.0*fvals[i+2][1] + 16.0*fvals[i+3][1] - 3.0*fvals[i+4][1]) / (12.0*h);
    else if(i >= 2 && i < fvals.size() - 2)//5-point (4th order) midpoint finite difference
    {
      fp = (fvals[i-2][1] - 8.0*fvals[i-1][1] + 8.0*fvals[i+1][1] - fvals[i+2][1]) / (12.0*h);
      fpp = (fvals[i-1][1] - 2.0*fvals[i][1] + fvals[i+1][1]) / (h*h);//3-point (2nd order) for 2nd derivative
    }
    else//5-point (4th order) right endpoint finite difference
      fp = (-25.0*fvals[i][1] + 48.0*fvals[i-1][1] - 36.0*fvals[i-2][1] + 16.0*fvals[i-3][1] - 3.0*fvals[i-4][1]) / (12.0*h);
    fvals[i][2] = fp;
    fvals[i][3] = fpp;
  }
  //Fill in endpoint values for 2nd derivatives
  fvals[0][3] = (-25.0*fvals[0][2] + 48.0*fvals[1][2] - 36.0*fvals[2][2] + 16.0*fvals[3][2] - 3.0*fvals[4][2]) / (12.0*h);
  fvals[1][3] = (-25.0*fvals[1][2] + 48.0*fvals[2][2] - 36.0*fvals[3][2] + 16.0*fvals[4][2] - 3.0*fvals[5][2]) / (12.0*h);
  int end = fvals.size() - 1;
  fvals[end-1][3] = (-25.0*fvals[end-1][2] + 48.0*fvals[end-2][2] - 36.0*fvals[end-3][2] + 16.0*fvals[end-4][2] - 3.0*fvals[end-5][2]) / (12.0*h);
  fvals[end][3] = (-25.0*fvals[end][2] + 48.0*fvals[end-1][2] - 36.0*fvals[end-2][2] + 16.0*fvals[i-3][2] - 3.0*fvals[i-4][2]) / (12.0*h);
  return fvals;
}

/*Function to transform vector of 4-tuples (z, f(z) = w, f'(z), f''(z)) in Kruskal coordinates to
a vector of 4-tuples (rho, g(rho) = tau, g'(rho), g''(rho)) in Tolman-Bondi coordinates in the
Glue and Friedman regions of the metric
*/
vector< array<double,4> > gGenerator(vector< array<double,4> > f)
{
  vector< array<double,4> > g;//Values of (rho, g(rho) = tau, g'(rho), g''(rho))
  int i;
  for(i = 0; i < f.size(); i++)
  {
    array<double,2> zw = {f[i][0], f[i][1]};
    array<double,2> rhotau = KTBTransform2(zw);
    if(rhotau[0] >= BETA)//Glue and Friedman regions of metric
    {
      array<double,4> gElement = {rhotau[0], rhotau[1], 0.0, 0.0};
      g.push_back(gElement);
    }
  }
  g = finiteDiff2(g);
  return g;
}

/*Function to calculate the mean curvature H in the Gluing and Friedman region from
  4-tuple (z, f(z) = w, f'(z), f''(z)) represented by kcoords*/
double calcHGF2(array<double,4> tbcoords)
{
  double r = tbcoords[0];
  double t = tbcoords[1];
  double fp = tbcoords[2];
  double fpp = tbcoords[3];
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
vector< array<double,2> > calculateH2(vector< array<double,4> > fVals)
{
  vector< array<double,2> > hValues;
  vector< array<double,4> > gVals = gGenerator(fVals);
  int i;
  int i_g_start = fVals.size() - gVals.size();//Where the gluing region of the metric starts
  double H;
  for(i = 0; i < fVals.size(); i++)
  {
    array<double,4> fElement = fVals[i];
    if(i < i_g_start)//Schwarzschild region, rho <= 0.75
      H = calcHS(fElement);
    else//Glue and Friedman Regions
    {
      array<double,4> gElement = gVals[i - i_g_start];
      H = calcHGF2(gElement);
    }
    array<double,2> nElem = {fElement[0], H};//Build ordered pair
    hValues.push_back(nElem);//Append ordered pair to vector
  }
  return hValues;
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
    double alpha = 1.0 / sqrt((32.0 / r)*exp(-r/2.0));
    double df_ds = (alpha/sqrt(1.0 - (fp_cur*fp_cur))) - sqrt(((alpha*alpha)/(1 - (fp_cur*fp_cur))) - 1.0);
    array<double,2> zw = {z_cur, f_cur};
    array<double,2> rt = KTBTransform2(zw);
    printf("rho: %lf\n", rt[0]);
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
  int numEvolutions = 100;
  vector< array<double,4> > f_test = evolve(numEvolutions);
  outputToFile4(f_test);
}

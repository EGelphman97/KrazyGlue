/*Eric Gelphman
  University of California, San Diego(UCSD)
  Department of Mathematics
  Irwin and Joan Jacobs School of Engineering Department of Electrical and Computer Engineering(ECE)

  December 31, 2017
  fulminePlusPlus, a C++ extension that performs all mathematical and file I/0 operations for KG2

  Version 1.0.1
*/

#ifdef __cplusplus
#include <Python.h>//Header file needed to interface with Python
#include <stdio.h>
#include <vector>//std::vector
#include <array>//std::array
#include <cmath>
using namespace std;

double ALPHA = 0.1;//rho parameters
double BETA = 0.5;
double GAMMA = 2.0;
double EPSILON = 1.0e-15;//Slightly better than machine precision

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

/*Helper function to tranfrom an input 4-tuple (z, f(z) = w, f'(z), f''(z)) in Kruskal coordinates to
an output 4-tuple (rho, g(rho) = tau, g'(rho), g''(rho)) in Tolman-Bondi Coordinates*/
array<double,4> KTBTransform(array<double,4> vals)
{
  array<double,4> result;
  double z = vals[0];
  double w = vals[1];//f(z)
  double d1f = vals[2];//f'(z)
  double d2f = vals[3];//f''(z)
  double U = w - z;//U = w - z
  double V = w + z;//V = w + z
  double t = 2.0*log(-V / U);
  double r = 2.0*LambertW((-U*V) / exp(1)) + 2.0;
  double tau = 2.0*sqrt(2.0)*sqrt(r) + t - 2.0*log(fabs((sqrt(2.0) + sqrt(r))/(sqrt(2.0) - sqrt(r))));
  double rho = tau + sqrt((2.0/9.0)*pow(r,3));
  result[0] = rho;
  result[1] = tau;//g(rho)
  double r1 = pow(4.5, 1.0/3.0)*pow(rho - tau, 2.0/3.0);
  double test = (sqrt(2.0) + sqrt(r1)) / (sqrt(2.0 - sqrt(r1)));//term in absolute value in expressions for g'(rho), g''(rho)
  if(test >= 0.0)//Positive domain of absolute value term
  {
    //g'(rho)
    result[2] = (-1.650963624447314*(exp(tau/2.)*(-1.9020315595498223e-16 + 0.27516060407455223*pow(rho,2) + 9.510157797749112e-17*pow(rho - tau,0.3333333333333333) - 0.6666666666666667*rho*pow(rho - tau,0.3333333333333333) + 0.40380457618491983*pow(rho - tau,0.6666666666666666) - 0.5503212081491045*rho*tau + 0.6666666666666667*pow(rho - tau,0.3333333333333333)*tau +
                0.27516060407455223*pow(tau,2) + d1f*(1.9020315595498223e-16 - 0.27516060407455223*pow(rho,2) - 9.510157797749112e-17*pow(rho - tau,0.3333333333333333) + 0.6666666666666667*rho*pow(rho - tau,0.3333333333333333) - 0.40380457618491983*pow(rho - tau,0.6666666666666666) + 0.5503212081491045*rho*tau - 0.6666666666666667*pow(rho - tau,0.3333333333333333)*tau -
                0.27516060407455223*pow(tau,2))) + exp(1.8171205928321397*pow(rho - tau,0.3333333333333333))*(-1.902031559549822e-16 + 0.2751606040745523*pow(rho,2) - 9.51015779774911e-17*pow(rho - tau,0.3333333333333333) + 0.4038045761849201*pow(rho - tau,0.6666666666666666) + rho*(9.510157797749112e-17 - 0.6666666666666667*pow(rho - tau,0.3333333333333333) +
                4.755078898874555e-17*pow(rho - tau,0.6666666666666666) - 0.5503212081491046*tau) - 9.51015779774911e-17*tau + 0.6666666666666667*pow(rho - tau,0.3333333333333333)*tau - 4.755078898874555e-17*pow(rho - tau,0.6666666666666666)*tau + 0.2751606040745523*pow(tau,2) + d1f*(-1.902031559549822e-16 + 0.2751606040745523*pow(rho,2) - 9.51015779774911e-17*pow(rho - tau,0.3333333333333333) +
                0.4038045761849201*pow(rho - tau,0.6666666666666666) + rho*(9.510157797749112e-17 - 0.6666666666666667*pow(rho - tau,0.3333333333333333) + 4.755078898874555e-17*pow(rho - tau,0.6666666666666666) - 0.5503212081491046*tau) + (-9.510157797749112e-17 + 0.6666666666666667*pow(rho - tau,0.3333333333333333) - 4.755078898874555e-17*pow(rho - tau,0.6666666666666666))*tau +
                0.2751606040745523*pow(tau,2)))))/(exp(1.8171205928321397*pow(rho - tau,0.3333333333333333))*(-0.4127409061118284*(1. + 1.*d1f)*pow(rho,2)*pow(rho - tau,0.3333333333333333) - 2.355138688025663e-16*pow(rho - tau,0.6666666666666666) + 0.6057068642773802*tau - 1.0000000000000004*pow(rho - tau,0.6666666666666666)*tau -
                0.4127409061118284*pow(rho - tau,0.3333333333333333)*pow(tau,2) + d1f*(-2.355138688025663e-16*pow(rho - tau,0.6666666666666666) + 0.6057068642773802*tau - 1.0000000000000004*pow(rho - tau,0.6666666666666666)*tau - 0.4127409061118284*pow(rho - tau,0.3333333333333333)*pow(tau,2)) + rho*(-0.6057068642773802 + 1.0000000000000004*pow(rho - tau,0.6666666666666666) +
                0.8254818122236568*pow(rho - tau,0.3333333333333333)*tau + d1f*(-0.6057068642773802 + 1.0000000000000004*pow(rho - tau,0.6666666666666666) + 0.8254818122236568*pow(rho - tau,0.3333333333333333)*tau))) + exp(tau/2.)*(-1.5700924586837754e-16 + pow(rho,2)*(1.5700924586837754e-16 + d1f*(-1.5700924586837754e-16 - 0.4127409061118284*pow(rho - tau,0.3333333333333333)) +
                0.4127409061118284*pow(rho - tau,0.3333333333333333)) + (-0.6057068642773802 + d1f*(0.6057068642773802 + 7.850462293418877e-17*pow(rho - tau,0.3333333333333333) - 1.*pow(rho - tau,0.6666666666666666)) - 7.850462293418877e-17*pow(rho - tau,0.3333333333333333) + 1.*pow(rho - tau,0.6666666666666666))*tau + (1.5700924586837754e-16 + d1f*(-1.5700924586837754e-16 -
                0.4127409061118284*pow(rho - tau,0.3333333333333333)) + 0.4127409061118284*pow(rho - tau,0.3333333333333333))*pow(tau,2) + rho*(0.6057068642773802 + 7.850462293418877e-17*pow(rho - tau,0.3333333333333333) - 1.*pow(rho - tau,0.6666666666666666) - 3.140184917367551e-16*tau - 0.8254818122236568*pow(rho - tau,0.3333333333333333)*tau + d1f*(-0.6057068642773802 -
                7.850462293418877e-17*pow(rho - tau,0.3333333333333333) + 1.*pow(rho - tau,0.6666666666666666) + 3.140184917367551e-16*tau + 0.8254818122236568*pow(rho - tau,0.3333333333333333)*tau))));
    //g''(rho)
    result[3] = (-1.650963624447314*d2f*pow(1.100642416298209 - 1.*pow(rho - tau,0.3333333333333333),2)*(1.100642416298209 + 1.*pow(rho - tau,0.3333333333333333))*(-1.2114137285547597 + 1.*pow(rho - tau,0.6666666666666666))*((0.3668808054327363*exp(tau/2.))/pow(1.100642416298209 -
                1.*pow(rho - tau,0.3333333333333333),2) - (0.3668808054327363*exp(1.8171205928321397*pow(rho - tau,0.3333333333333333)))/(1.2114137285547597 - 1.*pow(rho - tau,0.6666666666666666)) + (0.3028534321386899*exp(1.8171205928321397*pow(rho - tau,0.3333333333333333))*(-1.2114137285547597 + 0.9085602964160698*rho +
                1.*pow(rho - tau,0.6666666666666666) - 0.9085602964160698*tau))/(-1.2114137285547597 + 1.*pow(rho - tau,0.6666666666666666)) + (0.3028534321386899*exp(tau/2.)*(1.100642416298209 + 1.*pow(rho - tau,0.3333333333333333))*(-1.2114137285547597 - 0.9085602964160698*rho + 1.*pow(rho - tau,0.6666666666666666) +
                0.9085602964160698*tau))/((-1.100642416298209 + 1.*pow(rho - tau,0.3333333333333333))*(-1.2114137285547597 + 1.*pow(rho - tau,0.6666666666666666))) - (1.0000000000000004*(exp(1.8171205928321397*pow(rho - tau,0.3333333333333333))*(pow(rho,2)*(7.850462293418876e-17 - 0.41274090611182834*pow(rho - tau,0.3333333333333333)) -
                7.850462293418876e-17*pow(rho - tau,0.6666666666666666) + (0.6057068642773801 + 7.850462293418876e-17*pow(rho - tau,0.3333333333333333) - 1.*pow(rho - tau,0.6666666666666666))*tau + (7.850462293418876e-17 - 0.41274090611182834*pow(rho - tau,0.3333333333333333))*pow(tau,2) + rho*(-0.6057068642773801 - 7.850462293418876e-17*pow(rho - tau,0.3333333333333333) +
                1.*pow(rho - tau,0.6666666666666666) + (-1.5700924586837752e-16 + 0.8254818122236567*pow(rho - tau,0.3333333333333333))*tau)) + exp(tau/2.)*(1.5700924586837752e-16 + pow(rho,2)*(-7.850462293418876e-17 - 0.41274090611182834*pow(rho - tau,0.3333333333333333)) + (0.6057068642773801 - 7.850462293418876e-17*pow(rho - tau,0.3333333333333333) -
                1.*pow(rho - tau,0.6666666666666666))*tau + (-7.850462293418876e-17 - 0.41274090611182834*pow(rho - tau,0.3333333333333333))*pow(tau,2) + rho*(-0.6057068642773801 + 7.850462293418876e-17*pow(rho - tau,0.3333333333333333) + 1.*pow(rho - tau,0.6666666666666666) + (1.5700924586837752e-16 +
                0.8254818122236567*pow(rho - tau,0.3333333333333333))*tau)))*(exp(tau/2.)*(-1.9020315595498223e-16 + 0.27516060407455223*pow(rho,2) + 9.510157797749112e-17*pow(rho - tau,0.3333333333333333) - 0.6666666666666667*rho*pow(rho - tau,0.3333333333333333) + 0.40380457618491983*pow(rho - tau,0.6666666666666666) - 0.5503212081491045*rho*tau +
                0.6666666666666667*pow(rho - tau,0.3333333333333333)*tau + 0.27516060407455223*pow(tau,2) + d1f*(1.9020315595498223e-16 - 0.27516060407455223*pow(rho,2) - 9.510157797749112e-17*pow(rho - tau,0.3333333333333333) + 0.6666666666666667*rho*pow(rho - tau,0.3333333333333333) - 0.40380457618491983*pow(rho - tau,0.6666666666666666) + 0.5503212081491045*rho*tau -
                0.6666666666666667*pow(rho - tau,0.3333333333333333)*tau - 0.27516060407455223*pow(tau,2))) + exp(1.8171205928321397*pow(rho - tau,0.3333333333333333))*(-1.902031559549822e-16 + 0.2751606040745523*pow(rho,2) - 9.51015779774911e-17*pow(rho - tau,0.3333333333333333) + 0.4038045761849201*pow(rho - tau,0.6666666666666666) + rho*(9.510157797749112e-17 -
                0.6666666666666667*pow(rho - tau,0.3333333333333333) + 4.755078898874555e-17*pow(rho - tau,0.6666666666666666) - 0.5503212081491046*tau) - 9.51015779774911e-17*tau + 0.6666666666666667*pow(rho - tau,0.3333333333333333)*tau - 4.755078898874555e-17*pow(rho - tau,0.6666666666666666)*tau + 0.2751606040745523*pow(tau,2) + d1f*(-1.902031559549822e-16 + 0.2751606040745523*pow(rho,2) -
                9.51015779774911e-17*pow(rho - tau,0.3333333333333333) + 0.4038045761849201*pow(rho - tau,0.6666666666666666) + rho*(9.510157797749112e-17 - 0.6666666666666667*pow(rho - tau,0.3333333333333333) + 4.755078898874555e-17*pow(rho - tau,0.6666666666666666) - 0.5503212081491046*tau) + (-9.510157797749112e-17 + 0.6666666666666667*pow(rho - tau,0.3333333333333333) - 4.755078898874555e-17*pow(rho - tau,0.6666666666666666))*tau +
                0.2751606040745523*pow(tau,2)))))/(pow(1.100642416298209 - 1.*pow(rho - tau,0.3333333333333333),2)*(1.100642416298209 + 1.*pow(rho - tau,0.3333333333333333))*(-1.2114137285547597 +
                1.*pow(rho - tau,0.6666666666666666))*(exp(1.8171205928321397*pow(rho - tau,0.3333333333333333))*(-0.4127409061118284*(1. + 1.*d1f)*pow(rho,2)*pow(rho - tau,0.3333333333333333) - 2.355138688025663e-16*pow(rho - tau,0.6666666666666666) + 0.6057068642773802*tau -
                1.0000000000000004*pow(rho - tau,0.6666666666666666)*tau - 0.4127409061118284*pow(rho - tau,0.3333333333333333)*pow(tau,2) + d1f*(-2.355138688025663e-16*pow(rho - tau,0.6666666666666666) + 0.6057068642773802*tau - 1.0000000000000004*pow(rho - tau,0.6666666666666666)*tau - 0.4127409061118284*pow(rho - tau,0.3333333333333333)*pow(tau,2)) + rho*(-0.6057068642773802 +
                1.0000000000000004*pow(rho - tau,0.6666666666666666) + 0.8254818122236568*pow(rho - tau,0.3333333333333333)*tau + d1f*(-0.6057068642773802 + 1.0000000000000004*pow(rho - tau,0.6666666666666666) + 0.8254818122236568*pow(rho - tau,0.3333333333333333)*tau))) + exp(tau/2.)*(-1.5700924586837754e-16 + pow(rho,2)*(1.5700924586837754e-16 + d1f*(-1.5700924586837754e-16 -
                0.4127409061118284*pow(rho - tau,0.3333333333333333)) + 0.4127409061118284*pow(rho - tau,0.3333333333333333)) + (-0.6057068642773802 + d1f*(0.6057068642773802 + 7.850462293418877e-17*pow(rho - tau,0.3333333333333333) - 1.*pow(rho - tau,0.6666666666666666)) - 7.850462293418877e-17*pow(rho - tau,0.3333333333333333) + 1.*pow(rho - tau,0.6666666666666666))*tau +
                (1.5700924586837754e-16 + d1f*(-1.5700924586837754e-16 - 0.4127409061118284*pow(rho - tau,0.3333333333333333)) + 0.4127409061118284*pow(rho - tau,0.3333333333333333))*pow(tau,2) + rho*(0.6057068642773802 + 7.850462293418877e-17*pow(rho - tau,0.3333333333333333) - 1.*pow(rho - tau,0.6666666666666666) - 3.140184917367551e-16*tau - 0.8254818122236568*pow(rho - tau,0.3333333333333333)*tau +
                d1f*(-0.6057068642773802 - 7.850462293418877e-17*pow(rho - tau,0.3333333333333333) + 1.*pow(rho - tau,0.6666666666666666) + 3.140184917367551e-16*tau +
                0.8254818122236568*pow(rho - tau,0.3333333333333333)*tau)))))))/(exp(1.8171205928321397*pow(rho - tau,0.3333333333333333))*(-0.4127409061118284*(1. + 1.*d1f)*pow(rho,2)*pow(rho - tau,0.3333333333333333) - 2.355138688025663e-16*pow(rho - tau,0.6666666666666666) + 0.6057068642773802*tau -
                1.0000000000000004*pow(rho - tau,0.6666666666666666)*tau - 0.4127409061118284*pow(rho - tau,0.3333333333333333)*pow(tau,2) + d1f*(-2.355138688025663e-16*pow(rho - tau,0.6666666666666666) + 0.6057068642773802*tau - 1.0000000000000004*pow(rho - tau,0.6666666666666666)*tau - 0.4127409061118284*pow(rho - tau,0.3333333333333333)*pow(tau,2)) + rho*(-0.6057068642773802 +
                1.0000000000000004*pow(rho - tau,0.6666666666666666) + 0.8254818122236568*pow(rho - tau,0.3333333333333333)*tau + d1f*(-0.6057068642773802 + 1.0000000000000004*pow(rho - tau,0.6666666666666666) + 0.8254818122236568*pow(rho - tau,0.3333333333333333)*tau))) + exp(tau/2.)*(-1.5700924586837754e-16 + pow(rho,2)*(1.5700924586837754e-16 + d1f*(-1.5700924586837754e-16 -
                0.4127409061118284*pow(rho - tau,0.3333333333333333)) + 0.4127409061118284*pow(rho - tau,0.3333333333333333)) + (-0.6057068642773802 + d1f*(0.6057068642773802 + 7.850462293418877e-17*pow(rho - tau,0.3333333333333333) - 1.*pow(rho - tau,0.6666666666666666)) - 7.850462293418877e-17*pow(rho - tau,0.3333333333333333) + 1.*pow(rho - tau,0.6666666666666666))*tau +
                (1.5700924586837754e-16 + d1f*(-1.5700924586837754e-16 - 0.4127409061118284*pow(rho - tau,0.3333333333333333)) + 0.4127409061118284*pow(rho - tau,0.3333333333333333))*pow(tau,2) + rho*(0.6057068642773802 + 7.850462293418877e-17*pow(rho - tau,0.3333333333333333) - 1.*pow(rho - tau,0.6666666666666666) - 3.140184917367551e-16*tau - 0.8254818122236568*pow(rho - tau,0.3333333333333333)*tau +
                d1f*(-0.6057068642773802 - 7.850462293418877e-17*pow(rho - tau,0.3333333333333333) + 1.*pow(rho - tau,0.6666666666666666) + 3.140184917367551e-16*tau + 0.8254818122236568*pow(rho - tau,0.3333333333333333)*tau))));
  }
  else if(test < 0.0)//Negative domain of absolute value epxression
  {
    //g'(rho)
    result[2] = (1.650963624447314*(exp(tau/2.)*(1.9020315595498223e-16 - 0.27516060407455223*pow(rho,2) - 9.510157797749112e-17*pow(rho - tau,0.3333333333333333) + 0.6666666666666667*rho*pow(rho - tau,0.3333333333333333) - 0.40380457618491983*pow(rho - tau,0.6666666666666666) + 0.5503212081491045*rho*tau -
                0.6666666666666667*pow(rho - tau,0.3333333333333333)*tau - 0.27516060407455223*pow(tau,2) + d1f*(-1.9020315595498223e-16 + 0.27516060407455223*pow(rho,2) + 9.510157797749112e-17*pow(rho - tau,0.3333333333333333) - 0.6666666666666667*rho*pow(rho - tau,0.3333333333333333) + 0.40380457618491983*pow(rho - tau,0.6666666666666666) -
                0.5503212081491045*rho*tau + 0.6666666666666667*pow(rho - tau,0.3333333333333333)*tau + 0.27516060407455223*pow(tau,2))) + exp(1.8171205928321397*pow(rho - tau,0.3333333333333333))*(-1.902031559549822e-16 + 0.2751606040745523*pow(rho,2) - 9.51015779774911e-17*pow(rho - tau,0.3333333333333333) + 0.4038045761849201*pow(rho - tau,0.6666666666666666) +
                rho*(9.510157797749112e-17 - 0.6666666666666667*pow(rho - tau,0.3333333333333333) + 4.755078898874555e-17*pow(rho - tau,0.6666666666666666) - 0.5503212081491046*tau) - 9.51015779774911e-17*tau + 0.6666666666666667*pow(rho - tau,0.3333333333333333)*tau - 4.755078898874555e-17*pow(rho - tau,0.6666666666666666)*tau + 0.2751606040745523*pow(tau,2) + d1f*(-1.902031559549822e-16 +
                0.2751606040745523*pow(rho,2) - 9.51015779774911e-17*pow(rho - tau,0.3333333333333333) + 0.4038045761849201*pow(rho - tau,0.6666666666666666) + rho*(9.510157797749112e-17 - 0.6666666666666667*pow(rho - tau,0.3333333333333333) + 4.755078898874555e-17*pow(rho - tau,0.6666666666666666) - 0.5503212081491046*tau) + (-9.510157797749112e-17 + 0.6666666666666667*pow(rho - tau,0.3333333333333333) -
                4.755078898874555e-17*pow(rho - tau,0.6666666666666666))*tau + 0.2751606040745523*pow(tau,2)))))/(exp(1.8171205928321397*pow(rho - tau,0.3333333333333333))*(0.4127409061118284*(1. + 1.*d1f)*pow(rho,2)*pow(rho - tau,0.3333333333333333) + 2.355138688025663e-16*pow(rho - tau,0.6666666666666666) - 0.6057068642773802*tau +
                1.0000000000000004*pow(rho - tau,0.6666666666666666)*tau + 0.4127409061118284*pow(rho - tau,0.3333333333333333)*pow(tau,2) + d1f*(2.355138688025663e-16*pow(rho - tau,0.6666666666666666) - 0.6057068642773802*tau + 1.0000000000000004*pow(rho - tau,0.6666666666666666)*tau + 0.4127409061118284*pow(rho - tau,0.3333333333333333)*pow(tau,2)) + rho*(0.6057068642773802 -
                1.0000000000000004*pow(rho - tau,0.6666666666666666) - 0.8254818122236568*pow(rho - tau,0.3333333333333333)*tau + d1f*(0.6057068642773802 - 1.0000000000000004*pow(rho - tau,0.6666666666666666) - 0.8254818122236568*pow(rho - tau,0.3333333333333333)*tau))) + exp(tau/2.)*(-1.5700924586837754e-16 + pow(rho,2)*(1.5700924586837754e-16 + d1f*(-1.5700924586837754e-16 -
                0.4127409061118284*pow(rho - tau,0.3333333333333333)) + 0.4127409061118284*pow(rho - tau,0.3333333333333333)) + (-0.6057068642773802 + d1f*(0.6057068642773802 + 7.850462293418877e-17*pow(rho - tau,0.3333333333333333) - 1.*pow(rho - tau,0.6666666666666666)) - 7.850462293418877e-17*pow(rho - tau,0.3333333333333333) + 1.*pow(rho - tau,0.6666666666666666))*tau +
                (1.5700924586837754e-16 + d1f*(-1.5700924586837754e-16 - 0.4127409061118284*pow(rho - tau,0.3333333333333333)) + 0.4127409061118284*pow(rho - tau,0.3333333333333333))*pow(tau,2) + rho*(0.6057068642773802 + 7.850462293418877e-17*pow(rho - tau,0.3333333333333333) - 1.*pow(rho - tau,0.6666666666666666) - 3.140184917367551e-16*tau - 0.8254818122236568*pow(rho - tau,0.3333333333333333)*tau +
                d1f*(-0.6057068642773802 - 7.850462293418877e-17*pow(rho - tau,0.3333333333333333) + 1.*pow(rho - tau,0.6666666666666666) + 3.140184917367551e-16*tau + 0.8254818122236568*pow(rho - tau,0.3333333333333333)*tau))));
    //g''(rho)
    result[3] = (1.650963624447314*d2f*pow(1.100642416298209 - 1.*pow(rho - tau,0.3333333333333333),2)*(1.100642416298209 + 1.*pow(rho - tau,0.3333333333333333))*(-1.2114137285547597 + 1.*pow(rho - tau,0.6666666666666666))*((-0.3668808054327363*exp(tau/2.))/pow(1.100642416298209 -
                1.*pow(rho - tau,0.3333333333333333),2) - (0.3668808054327363*exp(1.8171205928321397*pow(rho - tau,0.3333333333333333)))/(1.2114137285547597 - 1.*pow(rho - tau,0.6666666666666666)) + (0.3028534321386899*exp(1.8171205928321397*pow(rho - tau,0.3333333333333333))*(-1.2114137285547597 + 0.9085602964160698*rho + 1.*pow(rho - tau,0.6666666666666666) -
                0.9085602964160698*tau))/(-1.2114137285547597 + 1.*pow(rho - tau,0.6666666666666666)) - (0.3028534321386899*exp(tau/2.)*(1.100642416298209 + 1.*pow(rho - tau,0.3333333333333333))*(-1.2114137285547597 - 0.9085602964160698*rho + 1.*pow(rho - tau,0.6666666666666666) + 0.9085602964160698*tau))/((-1.100642416298209 +
                1.*pow(rho - tau,0.3333333333333333))*(-1.2114137285547597 + 1.*pow(rho - tau,0.6666666666666666))) - (1.0000000000000004*(exp(1.8171205928321397*pow(rho - tau,0.3333333333333333))*(pow(rho,2)*(-7.850462293418876e-17 + 0.41274090611182834*pow(rho - tau,0.3333333333333333)) + 7.850462293418876e-17*pow(rho - tau,0.6666666666666666) +
                (-0.6057068642773801 - 7.850462293418876e-17*pow(rho - tau,0.3333333333333333) + 1.*pow(rho - tau,0.6666666666666666))*tau + (-7.850462293418876e-17 + 0.41274090611182834*pow(rho - tau,0.3333333333333333))*pow(tau,2) + rho*(0.6057068642773801 + 7.850462293418876e-17*pow(rho - tau,0.3333333333333333) - 1.*pow(rho - tau,0.6666666666666666) +
                (1.5700924586837752e-16 - 0.8254818122236567*pow(rho - tau,0.3333333333333333))*tau)) + exp(tau/2.)*(1.5700924586837752e-16 + pow(rho,2)*(-7.850462293418876e-17 - 0.41274090611182834*pow(rho - tau,0.3333333333333333)) + (0.6057068642773801 - 7.850462293418876e-17*pow(rho - tau,0.3333333333333333) - 1.*pow(rho - tau,0.6666666666666666))*tau +
                (-7.850462293418876e-17 - 0.41274090611182834*pow(rho - tau,0.3333333333333333))*pow(tau,2) + rho*(-0.6057068642773801 + 7.850462293418876e-17*pow(rho - tau,0.3333333333333333) + 1.*pow(rho - tau,0.6666666666666666) + (1.5700924586837752e-16 + 0.8254818122236567*pow(rho - tau,0.3333333333333333))*tau)))*(exp(tau/2.)*(1.9020315595498223e-16 -
                0.27516060407455223*pow(rho,2) - 9.510157797749112e-17*pow(rho - tau,0.3333333333333333) + 0.6666666666666667*rho*pow(rho - tau,0.3333333333333333) - 0.40380457618491983*pow(rho - tau,0.6666666666666666) + 0.5503212081491045*rho*tau - 0.6666666666666667*pow(rho - tau,0.3333333333333333)*tau - 0.27516060407455223*pow(tau,2) + d1f*(-1.9020315595498223e-16 +
                0.27516060407455223*pow(rho,2) + 9.510157797749112e-17*pow(rho - tau,0.3333333333333333) - 0.6666666666666667*rho*pow(rho - tau,0.3333333333333333) + 0.40380457618491983*pow(rho - tau,0.6666666666666666) - 0.5503212081491045*rho*tau + 0.6666666666666667*pow(rho - tau,0.3333333333333333)*tau + 0.27516060407455223*pow(tau,2))) +
                exp(1.8171205928321397*pow(rho - tau,0.3333333333333333))*(-1.902031559549822e-16 + 0.2751606040745523*pow(rho,2) - 9.51015779774911e-17*pow(rho - tau,0.3333333333333333) + 0.4038045761849201*pow(rho - tau,0.6666666666666666) + rho*(9.510157797749112e-17 - 0.6666666666666667*pow(rho - tau,0.3333333333333333) + 4.755078898874555e-17*pow(rho - tau,0.6666666666666666) - 0.5503212081491046*tau) -
                9.51015779774911e-17*tau + 0.6666666666666667*pow(rho - tau,0.3333333333333333)*tau - 4.755078898874555e-17*pow(rho - tau,0.6666666666666666)*tau + 0.2751606040745523*pow(tau,2) + d1f*(-1.902031559549822e-16 + 0.2751606040745523*pow(rho,2) - 9.51015779774911e-17*pow(rho - tau,0.3333333333333333) + 0.4038045761849201*pow(rho - tau,0.6666666666666666) + rho*(9.510157797749112e-17 -
                0.6666666666666667*pow(rho - tau,0.3333333333333333) + 4.755078898874555e-17*pow(rho - tau,0.6666666666666666) - 0.5503212081491046*tau) + (-9.510157797749112e-17 + 0.6666666666666667*pow(rho - tau,0.3333333333333333) - 4.755078898874555e-17*pow(rho - tau,0.6666666666666666))*tau + 0.2751606040745523*pow(tau,2)))))/(pow(1.100642416298209 -
                1.*pow(rho - tau,0.3333333333333333),2)*(1.100642416298209 + 1.*pow(rho - tau,0.3333333333333333))*(-1.2114137285547597 + 1.*pow(rho - tau,0.6666666666666666))*(exp(1.8171205928321397*pow(rho - tau,0.3333333333333333))*(0.4127409061118284*(1. + 1.*d1f)*pow(rho,2)*pow(rho - tau,0.3333333333333333) +
                2.355138688025663e-16*pow(rho - tau,0.6666666666666666) - 0.6057068642773802*tau + 1.0000000000000004*pow(rho - tau,0.6666666666666666)*tau + 0.4127409061118284*pow(rho - tau,0.3333333333333333)*pow(tau,2) + d1f*(2.355138688025663e-16*pow(rho - tau,0.6666666666666666) - 0.6057068642773802*tau + 1.0000000000000004*pow(rho - tau,0.6666666666666666)*tau +
                0.4127409061118284*pow(rho - tau,0.3333333333333333)*pow(tau,2)) + rho*(0.6057068642773802 - 1.0000000000000004*pow(rho - tau,0.6666666666666666) - 0.8254818122236568*pow(rho - tau,0.3333333333333333)*tau + d1f*(0.6057068642773802 - 1.0000000000000004*pow(rho - tau,0.6666666666666666) - 0.8254818122236568*pow(rho - tau,0.3333333333333333)*tau))) +
                exp(tau/2.)*(-1.5700924586837754e-16 + pow(rho,2)*(1.5700924586837754e-16 + d1f*(-1.5700924586837754e-16 - 0.4127409061118284*pow(rho - tau,0.3333333333333333)) + 0.4127409061118284*pow(rho - tau,0.3333333333333333)) + (-0.6057068642773802 + d1f*(0.6057068642773802 + 7.850462293418877e-17*pow(rho - tau,0.3333333333333333) -
                1.*pow(rho - tau,0.6666666666666666)) - 7.850462293418877e-17*pow(rho - tau,0.3333333333333333) + 1.*pow(rho - tau,0.6666666666666666))*tau + (1.5700924586837754e-16 + d1f*(-1.5700924586837754e-16 - 0.4127409061118284*pow(rho - tau,0.3333333333333333)) + 0.4127409061118284*pow(rho - tau,0.3333333333333333))*pow(tau,2) + rho*(0.6057068642773802 +
                7.850462293418877e-17*pow(rho - tau,0.3333333333333333) - 1.*pow(rho - tau,0.6666666666666666) - 3.140184917367551e-16*tau - 0.8254818122236568*pow(rho - tau,0.3333333333333333)*tau + d1f*(-0.6057068642773802 - 7.850462293418877e-17*pow(rho - tau,0.3333333333333333) + 1.*pow(rho - tau,0.6666666666666666) + 3.140184917367551e-16*tau +
                0.8254818122236568*pow(rho - tau,0.3333333333333333)*tau)))))))/(exp(1.8171205928321397*pow(rho - tau,0.3333333333333333))*(0.4127409061118284*(1. + 1.*d1f)*pow(rho,2)*pow(rho - tau,0.3333333333333333) + 2.355138688025663e-16*pow(rho - tau,0.6666666666666666) - 0.6057068642773802*tau + 1.0000000000000004*pow(rho - tau,0.6666666666666666)*tau +
                0.4127409061118284*pow(rho - tau,0.3333333333333333)*pow(tau,2) + d1f*(2.355138688025663e-16*pow(rho - tau,0.6666666666666666) - 0.6057068642773802*tau + 1.0000000000000004*pow(rho - tau,0.6666666666666666)*tau + 0.4127409061118284*pow(rho - tau,0.3333333333333333)*pow(tau,2)) + rho*(0.6057068642773802 - 1.0000000000000004*pow(rho - tau,0.6666666666666666) -
                0.8254818122236568*pow(rho - tau,0.3333333333333333)*tau + d1f*(0.6057068642773802 - 1.0000000000000004*pow(rho - tau,0.6666666666666666) - 0.8254818122236568*pow(rho - tau,0.3333333333333333)*tau))) + exp(tau/2.)*(-1.5700924586837754e-16 + pow(rho,2)*(1.5700924586837754e-16 + d1f*(-1.5700924586837754e-16 - 0.4127409061118284*pow(rho - tau,0.3333333333333333)) +
                0.4127409061118284*pow(rho - tau,0.3333333333333333)) + (-0.6057068642773802 + d1f*(0.6057068642773802 + 7.850462293418877e-17*pow(rho - tau,0.3333333333333333) - 1.*pow(rho - tau,0.6666666666666666)) - 7.850462293418877e-17*pow(rho - tau,0.3333333333333333) + 1.*pow(rho - tau,0.6666666666666666))*tau + (1.5700924586837754e-16 + d1f*(-1.5700924586837754e-16 -
                0.4127409061118284*pow(rho - tau,0.3333333333333333)) + 0.4127409061118284*pow(rho - tau,0.3333333333333333))*pow(tau,2) + rho*(0.6057068642773802 + 7.850462293418877e-17*pow(rho - tau,0.3333333333333333) - 1.*pow(rho - tau,0.6666666666666666) - 3.140184917367551e-16*tau - 0.8254818122236568*pow(rho - tau,0.3333333333333333)*tau + d1f*(-0.6057068642773802 -
                7.850462293418877e-17*pow(rho - tau,0.3333333333333333) + 1.*pow(rho - tau,0.6666666666666666) + 3.140184917367551e-16*tau + 0.8254818122236568*pow(rho - tau,0.3333333333333333)*tau))));
  }
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
  else
  {
    array<double,2> indices = binarySearch(vals, key);
    val = (vals[indices[0]][1] + vals[indices[1]][1]) / 2.0;
  }
  return val;//Return average of w-values at the two indices in vals
}

/*Function to generate a vector of ordered pairs(implemented as arrays of size 2) (z, w)
with an even z-step given by the value of the variable z_step. Variable nPoints
is needed for the initial generation of rho-values*/
vector< array<double,2> > getConstZStep(double z_step, int nPoints)
{
  vector< array<double,2> > result;
  //Generate values of rho and tau for ALPHA < rho <= GAMMA
  double rho_step = (GAMMA - ALPHA) / (double)nPoints;
  double rho = ALPHA;
  double tau;
  vector< array<double,2> > zw;
  int i;
  for(i = 0; i < nPoints; i++)//Iterate based on list index i, increment rho seperately
  {
    if(rho < BETA)//Glue region
      tau = -1.84115*pow(rho,3) + 0.720532*pow(rho,2) + 0.6651*rho - 1.713532;
    else//Friedman region
      tau = -(double)37.0 / 30.0;
    double r = pow(4.5, 1.0/3.0)*pow(rho - tau, 2.0/3.0);
    double sqr = sqrt(r);
    double a = r + 2.0*log((r/2.0) - 1);
    double e_t_4 = exp(-(sqr/sqrt(2.0)) + tau/4.0)*sqrt(fabs((sqrt(2.0) + sqr)/(sqrt(2.0) - sqr)));
    double U = (-exp(a/4.0)) / e_t_4;
    double V = exp(a/4.0) * e_t_4;
    array<double,2> element = {(V - U) / 2.0, (V + U) / 2.0};// (z, f(z) = w)
    zw.push_back(element);//Append array to vector
    rho += rho_step;
  }
  double z = zw[0][0];
  for(i = 0; i < nPoints; i++)//Iterate based on index in vector
  {
    array<double,2> element;
    if(i == 0)//First element is special case
    {
      element[0] = z;
      element[1] = zw[0][1];//Corresponding w value for z0 = zw[0][0]
      result.push_back(element);//z, element at index 0 in zw
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
  int j = result.size() - 1;
  while(fabs(result[j][1] - result[j - 1][1]) <= EPSILON)//Trim end
  {
      result.pop_back();
      j--;
  }
  return result;
}

/*Function to compute the transfromation T: (rho, g(rho) = tau, g'(rho), g''(rho)) |--> (z, f(z) = w, f'(z), f''(z))
Returns a vector of arrays of size 4, with the array of size 4 representing the 4-tuple (z, f(z) = w, f'(z), f''(z))*/
vector< array<double,4> > TBKTransform(double z_step, int nPoints)
{
  vector< array<double,2> > constZStep = getConstZStep(z_step, nPoints);
  vector< array<double,4> > result;
  int i;
  double h = (constZStep[constZStep.size() - 1][0] - constZStep[0][0]) / 2.0;
  for(i = 0; i < constZStep.size(); i++)
  {
    array<double,4> element;
    double z = constZStep[i][0];
    double w = constZStep[i][1];
    double fp, fpp;//f'(z), f''(z)
    if(i == 0)//First endpoint
    {
      fp = (-3.0*w + 4.0*constZStep[i + 1][1] - constZStep[i + 2][1]) / (2.0*h);//3 - point endpoint
      double fpplus1 = (constZStep[i + 2][1] - w) / (2.0*h);
      double fpplus2 = (constZStep[i + 3][1] - constZStep[i + 1][1]) / (2.0*h);
      fpp = (-3.0*fp + 4.0*fpplus1 - fpplus2) / (2.0*h);//Compute f'' in similar manner to f'
    }
    else if(i == constZStep.size() - 1)//Last endpoint
    {
      fp = (-3.0*w + 4.0*constZStep[i - 1][1] - constZStep[i - 2][1]) / (2.0*h);//3 - point endpoint
      double fpminus1 = (w - constZStep[i - 2][1]) / (2.0*h);
      double fpminus2 = (constZStep[i - 3][1] - constZStep[i - 1][1]) / (2.0*h);
      fpp = (-3.0*fp + 4.0*fpminus1 - fpminus2) / (2.0*h);//Compute f'' in similar manner to f'
    }
    else//Use 3-point midpoint formulas
    {
      fp = (constZStep[i + 1][1] - constZStep[i - 1][1]) / (2.0*h);
      fpp = (constZStep[i - 1][1] - 2.0*w + constZStep[i + 1][1]) / (h*h);
    }
    element[0] = z;
    element[1] = w;
    element[2] = fp;
    element[3] = fpp;
    result.push_back(element);
  }
  return result;
}

/*Function to convert a vector of stl arrays of size 4 to a single array,
with the first index storing the size of the array*/
double* vectorToArray(vector< array<double,4> > fVals)
{
  double arr[4*fVals.size() + 1];//First index stores
  int i,j,k;
  arr[0] = arr[4*fVals.size()];
  k = 1;
  for(i = 0; i < 4*fVals.size(); i++)
  {
    for(j = 0; j < 4; j++)
    {
      arr[k] = fVals[i][j];
      k++;
    }
  }
  return arr;
}

/*Function to generate an array of doubles representing values of z, w, and the first and second derivatives.
The firstEach set of 4 indices i, i+1, i+2, i+3represents the 4-tuple (z, f(z) = w, f'(z), f''(z)).*/
double* fGeneratorPlusPlus(double step_size, int nPoints)
{
  vector< array<double,4> > result;
  double a = 0.80360287983844092;//corresponding z-value for rho = ALPHA
  double z;
  for(z = 0.0; z <= a; z += step_size)
  {
      array<double, 4> element = {z, 0.0, 0.0, 0.0};
      result.push_back(element);
  }
  vector< array<double,4> > fGlueS = TBKTransform(step_size, nPoints);
  int i;
  for(i = 0; i < fGlueS.size(); i++)//Append contents of fGlueS to result
  {
    result.push_back(fGlueS[i]);
  }
  double* arr = vectorToArray(result);//Convert 4-tuples to array so it can be compiled using a C compiler
  return arr;//Return array
}

/*Function needed to expose the C++ function fGenerator()
Python parameters: float z_step representing step size for z,  int nPoints representing
number of rho values from which to generate z-values
Python return: array of tuples of size 4 of form (z, f(z) = w, f'(z), f''(z))
*/
static PyObject* fGenBuild(PyObject* self, PyObject* args)
{
  double z_step;
  int nPoints;
  if(!PyArg_ParseTuple(args, "di", &z_step, &nPoints))//If x is not a double, throw an exception
    return NULL;
  double* arr = fGeneratorPlusPlus(z_step, nPoints);//Perform numerical evaluation to obtain 4-tuple (z, f(z) = w, f'(z), f''(z))
  const unsigned tupleLength = 4;//Length of tuple, in terms of memory
  const unsigned arrayLength  = unsigned(arr[0]);//Length of array
  double* arr
  PyObject* list = PyListNew(0);//List of 4-tuples to be returned
  unsigned i;
  for(i = 0; i < arrayLength; i++)
  {
    PyObject* tuple = PyTupleNew(tupleLength);//Initialize tuple
    int j;
    for(j = 0; j < tupleLength; j++)
    {
      PyObject* val = PyFloatNew(arr[(4*i) + j])//Build array
      PyTuple_SET_ITEM(tuple, j, val)//Insert value into tuple
    }
    PyList_SET_ITEM(list, i, tuple)//Insert tuple into array
  }
  return list;
}

static char module_docstring[] =
    "fulminePlusPlus, a C++ module to perform numerical evaluation and differentiation for KG2";
static char evalFbarPara_docstring[] =
    "function to generate list of 4-tuples (z, f(z) = w, f'(z), f''(z)";

//Mapping of the names of Python methods to C functions
static PyMethodDef fulminePlusPlusMethods[] = {
      {"fGenerator5", fGenBuild, METH_VARARGS, fGenerator_docstring},//Mapping to fGeneratorPlusPlus()
      {NULL, NULL, 0, NULL}//Null terminator needed for formatting
};

//Module definition structure
static struct PyModuleDef fulminePlusPlusModule = {
    PyModuleDef_HEAD_INIT,
    "fulminePlusPlus",     //Name of module
    module_docstring,      //Documentation
    -1,
    fulminePlusPlusMethods
};

extern "C" {
#endif
//Function to initialize the module
PyMODINIT_FUNC PyInit_fulmine(void)
{
    return PyModule_Create(&fulminePlusPlusModule);
}
#ifdef __cplusplus
}
#endif
/*
//Function to calculate the mean curvature H(z,w) in the Schwarzschild region
double calcHS(double z, double f, double fp, double fpp)
{
  return 0.0;
}

//Function to calculate the mean curvature H(z,w) in the Gluing and Friedman region
double calcHGF(double z, double f, double fp, double fpp)
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
  vector< array<double,4> > fGen = fGeneratorPlusPlus(0.002, 1000);//fGenerator's Newest Form
  int i;
  for(i = 0; i < fGen.size(); i++)
  {
    array<double, 4> element = fGen[i];
    printf("{%lf,%lf}, ", element[0], element[1]);
  }
  return 0;
}
*/

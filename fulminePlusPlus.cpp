/*Eric Gelphman
  University of California, San Diego(UCSD)
  Department of Mathematics
  Irwin and Joan Jacobs School of Engineering Department of Electrical and Computer Engineering(ECE)

  January 1, 2018
  fulminePlusPlus, a C++ extension that performs all mathematical and file I/0 operations for KG2

  Version 1.1.0
*/

#include <stdio.h>
#include <iostream>
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
    array<double,2> rhotau = {rho, tau};// (rho, tau) to be transformed
    array<double,2> element = discTBKTransform(rhotau);//What is going to be returned
    zw.push_back(element);//Append array to vector
    rho += rho_step;
  }
  double z = zw[0][0];
  array<double,2> rhotauEnd = {GAMMA, -37.0 / 30.0};
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

/*Function to generate a vector of double arrays of size 4 representing values of z, w, and the first and second derivatives.
Each array represents the 4-tuple (z, f(z) = w, f'(z), f''(z)).*/
vector< array<double,4> > fGeneratorPlusPlus(double step_size, int nPoints)
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
  4-tuple (z, f(z) = w, f'(z), f''(z)) represented by kCoords*/
double calcHGF(array<double,4> kCoords)
{
  array<double,4> tbCoords = KTBTransform(kCoords);
  double r = tbCoords[0];
  double t = tbCoords[1];
  double fp = tbCoords[2];
  double fpp = tbCoords[3];
  double X, Y, dX_dr, dX_dt, dY_dr, dY_dt;
  static double H;
  if(r >= 1.2)//Friedman Region M(r) = r^3, t0(r) = 1, applies for all r >= 0.5
  {
    X = (pow(3,0.6666666666666666)*pow(r,2)*(1 - t))/(pow(2,0.3333333333333333)*pow(pow(r,6)*(1.0 - t),0.3333333333333333));
    dX_dr = 0.0;
    dX_dt = -((pow(2,0.6666666666666666)*pow(r,2))/(pow(3,0.3333333333333333)*pow(-(pow(r,6)*(-1.0 + t)),0.3333333333333333)));
    Y = pow(4.50, 1.0/3.0)*r*pow(1 - t, 2.0/3);
    dY_dr = pow(4.50, 1.0/3.0)*pow(1 - t, 2.0/3);
    dY_dt = (-1.100642416*r) / pow(1 - t, 1.0/3);
  }
  else//Glue Region applies for 0.1 < r < 0.5
  {
    X = (2*(-27 + 90*r - 93.75*pow(r,2) + 31.25*pow(r,3))*(1.512 - 21.6*r + 88.2*pow(r,2) - 138.25*pow(r,3) + 94.6875*pow(r,4) - 23.4375*pow(r,5)) +
        (-21.6 + 176.4*r - 414.75*pow(r,2) + 378.75*pow(r,3) - 117.1875*pow(r,4))*(6.4 - 27*r + 45*pow(r,2) - 31.25*pow(r,3) + 7.8125*pow(r,4) - t))/
        (pow(6,0.3333333333333333)*pow(pow(1.512 - 21.6*r + 88.2*pow(r,2) - 138.25*pow(r,3) + 94.6875*pow(r,4) -
        23.4375*pow(r,5),2)*(6.4 - 27*r + 45*pow(r,2) - 31.25*pow(r,3) + 7.8125*pow(r,4) - t),0.3333333333333333));
    double dX_dr_Num = 1.402168623996777 + 3.666382459143577e6*pow(r,18) - 684462.0742918561*pow(r,19) + 79767.01691360433*pow(r,20) - 4366.562320519018*pow(r,21) +
            pow(r,14)*(1.3939247798328805e8 - 278787.1283822919*t) + pow(r,16)*(3.829358379162845e7 - 10540.714840585912*t) - 0.19774845563650464*t - 0.008593825839646838*pow(t,2) +
            pow(r,17)*(-1.3742637198641371e7 + 739.4941234503592*t) + pow(r,15)*(-8.225734141495614e7 + 69305.14160522398*t) + pow(r,7)*(-1.3930824322660064e7 + 778518.5670809123*t -
            1382.8832281183295*pow(t,2)) + pow(r,9)*(-7.82507496515977e7 + 2.260651898757414e6*t - 1334.9178234153358*pow(t,2)) +
            pow(r,5)*(-939034.6828495631 + 86523.32391806085*t - 310.6557815252211*pow(t,2)) + pow(r,11)*(-1.8453572308964357e8 + 2.280388111738583e6*t - 270.5349708796563*pow(t,2)) +
            pow(r,3)*(-19222.543408078545 + 2497.581547038366*t - 8.752185508049237*pow(t,2)) + pow(r,13)*(-1.89129868163097e8 + 767132.6185375242*t - 5.503212081491041*pow(t,2)) +
            pow(r,2)*(1533.3627901175862 - 221.0799598868958*t - 0.059731696412250714*pow(t,2)) + r*(-70.92671275275359 + 10.645535348763083*t + 0.12584713841830317*pow(t,2)) +
            pow(r,12)*(2.073436365106003e8 - 1.529129577650688e6*t + 57.805739703981956*pow(t,2)) + pow(r,4)*(159294.30723341077 - 17788.78507195431*t + 73.94049325029711*pow(t,2)) +
            pow(r,10)*(1.3346535376715788e8 - 2.590495693345052e6*t + 744.2856542753781*pow(t,2)) + pow(r,6)*(4.1248685011178674e6 - 301727.5974644786*t + 804.0799110682765*pow(t,2)) +
            pow(r,8)*(3.696740152888469e7 - 1.5170766111478277e6*t + 1633.0799177990832*pow(t,2));
    double dX_dr_Denom = (pow(0.064512 - 0.9216000000000001*r + 3.7632000000000003*pow(r,2) - 5.898666666666666*pow(r,3) + 4.04*pow(r,4) - pow(r,5),2)*pow(pow(1.512 - 21.6*r + 88.2*pow(r,2) -
                         138.25*pow(r,3) + 94.6875*pow(r,4) - 23.4375*pow(r,5),2)*(6.4 - 27*r + 45*pow(r,2) - 31.25*pow(r,3) + 7.8125*pow(r,4) - t),0.3333333333333333)*(0.8192 -
                         3.456*r + 5.76*pow(r,2) - 4.0*pow(r,3) + pow(r,4) - 0.128*t));
    dX_dr = dX_dr_Num / dX_dr_Denom;
    dX_dt = (445.4140078080003 - 1.1174983978271486e8*pow(r,15) + 2.4200284481048584e7*pow(r,16) - 3.2347440719604488e6*pow(r,17) + 201165.67611694336*pow(r,18) + pow(r,8)*(1.2458457605343757e9 - 6.6827533453125e7*t) +
            pow(r,10)*(2.189734626064453e9 - 4.5997873388671875e7*t) + pow(r,6)*(2.6377120935000014e8 - 2.6682266362500004e7*t) + pow(r,12)*(1.4780344404907227e9 - 7.438005065917969e6*t) + pow(r,4)*(1.8044774044200003e7 - 2.756455191e6*t) +
            pow(r,14)*(3.5647639389038086e8 - 128746.03271484375*t) + pow(r,2)*(290505.19363200024 - 56618.245728000016*t) - 98.76142080000005*t + r*(-17265.847818240014 + 3628.3064832000014*t) +
            pow(r,3)*(-2.8311145168320006e6 + 497276.32752000005*t) + pow(r,13)*(-8.3333009765625e8 + 1.4563751220703125e6*t) + pow(r,5)*(-8.059053482640001e7 + 1.0264048002e7*t) + pow(r,11)*(-2.0313835853027346e9 + 2.2681193115234375e7*t) +
            pow(r,7)*(-6.524829870975003e8 + 4.960860266250001e7*t) + pow(r,9)*(-1.8609794385e9 + 6.5376476953125e7*t))/(3.0*pow(6,0.3333333333333333)*pow(pow(1.512 - 21.6*r + 88.2*pow(r,2) - 138.25*pow(r,3) + 94.6875*pow(r,4) -
            23.4375*pow(r,5),2)*(6.4 - 27*r + 45*pow(r,2) - 31.25*pow(r,3) + 7.8125*pow(r,4) - t),1.3333333333333333));
    Y = 1.6509636244473134*pow((1.512 - 21.6*r + 88.2*pow(r,2) - 138.25*pow(r,3) + 94.6875*pow(r,4) - 23.4375*pow(r,5))*pow(6.4 - 27*r + 45*pow(r,2) - 31.25*pow(r,3) +
        7.8125*pow(r,4) - t,2),0.3333333333333333);
    dY_dr = (-10234.130438716447*(0.8192 - 3.456*r + 5.76*pow(r,2) - 4.0*pow(r,3) + pow(r,4) - 0.128*t)*(0.09237551261538461 - 35.41070769230769*pow(r,5) + 21.409641025641026*pow(r,6) - 7.113846153846154*pow(r,7) + pow(r,8) +
            pow(r,2)*(7.27764676923077 - 0.17423753846153847*t) + pow(r,4)*(34.955421538461536 - 0.04923076923076923*t) + r*(-1.323625550769231 + 0.07410609230769231*t) + pow(r,3)*(-20.888024615384616 + 0.15911384615384616*t) -
            0.009074215384615385*t))/pow((1.512 - 21.6*r + 88.2*pow(r,2) - 138.25*pow(r,3) + 94.6875*pow(r,4) - 23.4375*pow(r,5))*pow(-6.4 + 27*r - 45*pow(r,2) + 31.25*pow(r,3) -
            7.8125*pow(r,4) + t,2),0.6666666666666666);
    dY_dt = (201.53364556241618*(-0.064512 + 0.9216000000000001*r - 3.7632000000000003*pow(r,2) + 5.898666666666666*pow(r,3) - 4.04*pow(r,4) + pow(r,5))*(0.8192 - 3.456*r + 5.76*pow(r,2) - 4.0*pow(r,3) + pow(r,4) - 0.128*t))/
            pow((1.512 - 21.6*r + 88.2*pow(r,2) - 138.25*pow(r,3) + 94.6875*pow(r,4) - 23.4375*pow(r,5))*pow(-6.4 + 27*r - 45*pow(r,2) + 31.25*pow(r,3) -
            7.8125*pow(r,4) + t,2),0.6666666666666666);
  }
  double H_Denom = X * Y * pow(pow(X,2) - pow(fp,2), 1.5);
  if(isnan(H_Denom))//If there is a negative number under the radical, throw an exception
  {
      printf("Error! Cannot Have a Negative Number Raised to the Power 3/2!\n");
      printf("r= %lf f'= %lf X= %lf\n", r, fp, X);
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
    if(z < 0.80360287983844092)//Schwarzschild region
    {
      H = calcHS(element);
    }
    else//Glue and Friedman Regions
    {
      H = calcHGF(element);
    }
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

int main()
{
  vector< array<double,4> > fGen = fGeneratorPlusPlus(0.002, 1000);//fGenerator's Newest Form
  outputToFileF(fGen);//Output 4-tuples (z, f(z) = w, f'(z), f''(z)) to file
  //vector< array<double,2> > h_vals = calculateH(fGen);
  //outputToFileH(h_vals);//Output ordered pairs (z, H) to file
  return 0;
}

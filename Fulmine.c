/*Eric Gelphman
  University of California San Diego(UCSD)
  Department of Mathematics
  Jacobs School of Engineering Department of Electrical and Computer Engineering(ECE)

  September 22, 2017

  Fulmine, a C module to perform arithmetic operations for fGenerator and KrazyGlue.
  Version 1.0.2
*/

#include <Python.h>//Header file needed to interface with Python
#include <stdio.h>
#include <math.h>

static char module_docstring[] =
    "Fulmine, a C module to perform numerical evaluation and differentiation for fGenerator and KrazyGlue";
static char evalFbarPara_docstring[] =
    "Function to evaluate r and fbar in the parametric region";
static char evalFbarExplicit_docstring[] =
    "Function to evaluate fbar in the glue and linear regions, where it is explicitly defined";
static char derivPara_docstring[] =
    "Function to calculate the first and second derivatives of fbar in the parametric region";
static char derivExplicit_docstring[] =
    "Function to calculate the first and second derivatives of fbar in the glue and linear regions, where it is explicitly defined";
static char calcHS_docstring[] = "Function to calculate the mean curvature H in the Schwarzschild region ";
static char calcHGF_docstring[] = "Function to calculate the mean curvature H in the Friedman and glue regions";

/*Function to evaluate r and fbar in the parametric region(r <= -1.0), where both
r and fbar(r) = t are defined in terms of the parameter x

Parameter x is the value of x that r and fbar are to evaluated at

Returns a pointer to an array which stores the values of r and fbar
*/
double* paraEvalFbar(double x)
{
    double a,r,fbar;
    double* result = malloc(sizeof(double)*2);//Manually allocate memory for array
    a = 1.0 - x;
    r = 4.0/(3.0*a) + 4.0*a - 2.0*log(2.0 - x) + 2.0*log(x);//Evaluation of r in the parametric region
    result[0] = r;
    if(r < -10.0)//6th order Taylor series approximation in parametric region for large negative r
    {
      fbar = -2.0*x - (5*pow(x,2))/2.0 - (35*pow(x,3))/12.0 - (105*pow(x,4))/32.0 - (231*pow(x,5))/64.0 - (1001*pow(x,6))/256.0;/* -
             (2145*pow(x,7))/512.0 - (36465*pow(x,8))/8192.0 - (230945*pow(x,9))/49152.0 - (323323*pow(x,10))/65536.0;*/
    }
    else//Regular evaluation for fbar
      fbar = ((double)4.0)/3.0 - 4.0/(3.0*pow(a,1.5));
    result[1] = fbar;
    return result;
}

/*Function to evaluate fbar in the glue(-1.0 < r < -0.8) and linear(r >= -0.8) regions,
where it is explicitly defined in terms of the independent variable r

Parameter r is the value of r that fbar is to be evaluated at

Returns a double which is the value of fbar for the input value of r
*/
double explicitEvalFbar(double r)
{
    double result = 0.0;
    if(r > 0.1 && r < 0.5)//Evaluation in glue region
      result = -1.84115*pow(r,3) + 0.720532*pow(r,2) - 0.33949*r - 0.380199;
    else if(r >= 0.5)//Evaluation in linear region
      result = -r - 0.1;
    return result;
}

/*Function to calculate the first and second order derivatives of fbar in
the parametric region

For r <= -50.0, taylor series of derivative is used for the numerical evaluation,
otherwise the derivative is evaluated using its explicit formula

Parameter x is the value of x the formula is to be evaluated at, r is the corresponding r-value

Returns a pointer to an array which stores the values of the first and second derivatives
*/
double* evalDerivativesPara(double x, double r)
{
    double* deriv = malloc(sizeof(double)*2);
    double d1, d2;//fbar', fbar''
    d1 = (x*(x -2))/(2*pow(1 - x, 4.5) - x*(x - 2));
    d2 = -((-2 + x)*pow(-1 + x,6)*x*(-4 - 10*x + 5*pow(x,2)))/(2.0*pow(2*sqrt(1 - x) + (2 - 8*sqrt(1 - x))*x +
          (-1 + 12*sqrt(1 - x))*pow(x,2) - 8*sqrt(1 - x)*pow(x,3) + 2*sqrt(1 - x)*pow(x,4),3));
    deriv[0] = d1;
    deriv[1] = d2;
    return deriv;
}

/*Function to calculate the first and second order derivatives of fbar in the
glue and linear regions, where fbar is defined expliticly. Therefore, the derivatives
are calculated explicitly using the power rule

Parameter r is the value of r the derivatives are to be evaluated at

Returns a pointer to an array which stores the values of the first and second derivatives
*/
double* evalDerivativesExplicit(double r)
{
    double dfbar_1, dfbar_2;//First, and second derivatives, respectively
    double* deriv = malloc(sizeof(double)*2);
    if(r > 0.1 && r < 0.5)//Glue region
    {
      dfbar_1 = -5.52345*pow(r,2) + 1.44106*r - 0.33949;
      dfbar_2 = -11.0469*r + 1.44106;
    }
    else//Linear region (r >= 0.5)
    {
      dfbar_1 = -1.0;
      dfbar_2 = 0.0;
    }
    deriv[0] = dfbar_1;
    deriv[1] = dfbar_2;
    return deriv;
}

/*Function to calculate H in the Schwarzschild region M(r) = 1, t0(r) = r, applies for all r <= 0.8

Paramters: r, fbar, fbarp, fpp are the values of r, fbar, fbar', and fbar'' from fGenerator

Return: Pointer to a double that stores the value of the mean curvature H
*/
double* calculateHS(double r, double fbar, double fbarp, double fpp)
{
  double X, Y, z, H_Num, H_Denom;
  static double H;
  z = -fbar + 4.0/3.0;//z = r - t, t = fbar + r - 4/3
  //printf("r: %.3lf fbar: %.60lf z: %.60lf\n", r, fbar, z);
  Y = pow(4.50, ((double)1.0)/3.0)*pow(z, ((double)2.0)/3.0);
  if(r <= -35.0)
  {
    double Xbar = ((fbar)/4.0 - (3*pow(-fbar,2))/8.0 - (53*pow(-fbar,3))/288.0 - (25*pow(-fbar,4))/768.0 - (41*pow(-fbar,5))/27648.0)/
                  (1.0 + 2.0*(-fbar) + (13*pow(-fbar,2))/9.0 + (65*pow(-fbar,3))/144.0 + (65*pow(-fbar,4))/1152.0 +
                  (13*pow(-fbar,5))/6912.0);
    X = Xbar + 1.0;
    H_Denom = X * Y * (pow(pow(Xbar,2) + 2.0*Xbar - pow(fbarp,2) - 2.0*fbarp, 1.5));
    double aleph = (2 + (7*(-fbar))/2.0 + (77*pow(-fbar,2))/36.0 + (77*pow(-fbar,3))/144.0 + (55*pow(-fbar,4))/1152.0 + (11*pow(-fbar,5))/13824.0)/
                   (1 + 2*(-fbar) + (13*pow(-fbar,2))/9.0 + (65*pow(-fbar,3))/144.0 + (65*pow(-fbar,4))/1152.0 + (13*pow(-fbar,5))/6912.);
    double alpha = (5 + (1609*(-fbar))/168.0 + (38305*pow(-fbar,2))/6048.0 + (27353*pow(-fbar,3))/16128.0 + (31081*pow(-fbar,4))/193536.0 + (407*pow(-fbar,5))/145152.0)/
                   (1 + (347*(-fbar))/168.0 + (2327*pow(-fbar,2))/1512.0 + (2665*pow(-fbar,3))/5376.0 + (6175*pow(-fbar,4))/96768.0 + (5083*pow(-fbar,5))/2.322432e6);
    double beta = (1.5 + (4647*(-fbar))/938.0 + (48949*pow(-fbar,2))/11256.0 + (499361*pow(-fbar,3))/360192.0 + (211511*pow(-fbar,4))/1.440768e6 + (5885*pow(-fbar,5))/2.161152e6)/
                  (1 + (3851*(-fbar))/1876.0 + (25675*pow(-fbar,2))/16884.0 + (263185*pow(-fbar,3))/540288.0 + (67405*pow(-fbar,4))/1.080576e6 + (55211*pow(-fbar,5))/2.5933824e7);
    double combined = (-fbar)/8.0 + (21*pow(-fbar,2))/64.0 - (437*pow(-fbar,3))/768.0 + (675*pow(-fbar,4))/1024.0 - (8119*pow(-fbar,5))/12288.0 + (30223*pow(-fbar,6))/49152.0 -
                      (322703*pow(-fbar,7))/589824.0 + (92891*pow(-fbar,8))/196608.0 - (16961825*pow(-fbar,9))/4.2467328e7 + (4181627*pow(-fbar,10))/1.2582912e7;
    H_Num = aleph*pow(fbarp, 3) + alpha*pow(fbarp, 2) + beta*fbarp + combined - 2.0*fpp;
    printf("r: %.3lf fbar' term: %.40lf constant: %.40lf f'' term:%.40lf\n",  r, beta*fbarp, combined, -2.0*fpp);
    printf("r: %.3lf Hd: %.40lf\n", r, H_Denom);
  }
  else
  {
    X = pow(2,((double)2.0)/3.0)/(pow(3,((double)1.0)/3.0)*pow(z,((double)1.0)/3.0));
    double dX_dr = -pow(2,((double)2.0)/3.0)/(3.0*pow(3,((double)1.0)/3.0)*pow(z,((double)4.0)/3.0));
    double fp = fbarp + 1.0;
    H_Denom = X * Y * pow(pow(X,2) - pow(fp,2), 1.5);
    double dX_dt = -dX_dr;
    double dY_dr = X;
    double dY_dt = -X;
    H_Num = (2.0*pow(fp,3)*dY_dr) - (pow(X,3)*Y*dX_dt) + ((X*Y*fp)*(dX_dr + 2.0*fp*dX_dt)) - (2.0*pow(X,4)*dY_dt) - (pow(X,2)*((Y*fpp)+(2.0*fp)*(dY_dr - fp*dY_dt)));
  }
  if(isnan(H_Denom))//If a negative number is under the radical, throw an exception
  {
    printf("Error! Cannot Have a Negative Number Raised to the Power 3/2!\n");
    printf("r= %lf X= %lf\n", r, X);
    return NULL;
  }
  else if(H_Denom == 0.0)//If the denominator = 0, throw an exception
  {
    printf("Error! Cannot Divide by 0!\n");
    printf("r= %lf\n", r);
    return NULL;
  }
  H = H_Num / H_Denom;
  printf("r: %.5lf fbar': %.15lf fbar'': %.15lf\n", r, fbarp, fpp);
  double* H_Ptr = &H;
  return H_Ptr;
}

/*
dX_dr = -0.25 - fbar/4.0 - (7*pow(-fbar,2))/32.0 + (35*pow(-fbar,3))/192.0 - (455*pow(-fbar,4))/3072.0 + (91*pow(-fbar,5))/768.0 -
        (1729*pow(-fbar,6))/18432.0 + (2717*pow(-fbar,7))/36864.0 - (67925*pow(-fbar,8))/1.179648e6 + (475475*pow(-fbar,9))/1.0616832e7 - (2947945*pow(-fbar,10))/8.4934656e7;
Xbar = fbar/4.0 + pow(-fbar,2)/8.0 - (7*pow(-fbar,3))/96.0 + (35*pow(-fbar,4))/768.0 - (91*pow(-fbar,5))/3072.0 +
        (91*pow(-fbar,6))/4608.0 - (247*pow(-fbar,7))/18432.0 + (2717*pow(-fbar,8))/294912.0 - (67925*pow(-fbar,9))/((double)1.0616832e7) + (95095*pow(-fbar,10))/((double)2.1233664e7);
alpha = 5.0 - (3.0*(-fbar))/4.0 + (3*pow(-fbar,2))/16.0 - pow(-fbar,3)/64.0 - (11*pow(-fbar,4))/256.0 + (61*pow(-fbar,5))/1024.0 -
                       (731*pow(-fbar,6))/12288.0 + (2609*pow(-fbar,7))/49152.0 - (8815*pow(-fbar,8))/196608.0 + (259741*pow(-fbar,9))/7.077888e6 - (833563*pow(-fbar,10))/2.8311552e7;
beta = 1.5 + (15*(-fbar))/8.0 - (57*pow(-fbar,2))/32.0 + (187*pow(-fbar,3))/128.0 - (589*pow(-fbar,4))/512.0 + (1823*pow(-fbar,5))/2048.0 -
                                     (16771*pow(-fbar,6))/24576.0 + (51145*pow(-fbar,7))/98304.0 - (155411*pow(-fbar,8))/393216.0 + (4239569*pow(-fbar,9))/1.4155776e7 -
                                     (12827387*pow(-fbar,10))/5.6623104e7;
*/

/* Function to calculate the numerical values of X, Y in the Glue and Friedman regions

Paramters: r, t, fp, fpp are the values of r, t, f'(r), and f''(r) from fGenerator(fbar converted to f)

Return: Pointer to a double that stores the value of the mean curvature H
*/
double* calculateHGF(double r, double t, double fp, double fpp)
{
    double X, Y, dX_dr, dX_dt, dY_dr, dY_dt;
    static double H;
    if(r >= 1.2)//Friedman Region M(r) = r^3, t0(r) = 1, applies for all r >= 1.2
    {
      X = (pow(3,0.6666666666666666)*pow(r,2)*(1 - t))/(pow(2,0.3333333333333333)*pow(pow(r,6)*(1.0 - t),0.3333333333333333));
      dX_dr = 0.0;
      dX_dt = -((pow(2,0.6666666666666666)*pow(r,2))/(pow(3,0.3333333333333333)*pow(-(pow(r,6)*(-1.0 + t)),0.3333333333333333)));
      Y = pow(4.50, 1.0/3.0)*r*pow(1 - t, 2.0/3);
      dY_dr = pow(4.50, 1.0/3.0)*pow(1 - t, 2.0/3);
      dY_dt = (-1.100642416*r) / pow(1 - t, 1.0/3);
    }
    else//Glue Region M(r) = 4.25r^3 - 7.35r^2 + 3.6r + 0.648, t0(r) = -1.25r^2 + 3r - 0.8, applies for 0.8 < r < 1.2
    {
      X = (2*(3 - 2.5*r)*(0.648 + 3.6*r - 7.35*pow(r,2) + 4.25*pow(r,3)) +
          (3.6 - 14.7*r + 12.75*pow(r,2))*(-0.8 + 3.0*r - 1.25*pow(r,2) - t))/
          (pow(6,0.3333333333333333)*pow(pow(0.648 + 3.6*r - 7.35*pow(r,2) +
          4.25*pow(r,3),2)*(-0.8 + 3.0*r - 1.25*pow(r,2) - t),0.3333333333333333));
      dX_dr = (0.3147944506877429 + 221.64186658205165*pow(r,10) - 27.28675990405975*pow(r,11) +
              pow(r,9)*(-768.8853224449925 - 28.066381615604318*t) + 0.573283980986305*t +
              0.28691491788233897*pow(t,2) + pow(r,8)*(1491.9154348970615 + 157.28180128901397*t) +
              pow(r,2)*(-11.748431131370616 - 23.240019899778524*t - 6.166611143524113*pow(t,2)) +
              pow(r,4)*(-12.129112454495935 + 56.34030249795842*t - 4.509139271939641*pow(t,2)) +
              pow(r,6)*(1207.2353490473233 + 481.2279501109907*t - 1.2880880741155186e-14*pow(t,2)) +
              pow(r,7)*(-1746.1824316337115 - 377.4219384095341*t - 1.6101100926443983e-15*pow(t,2)) +
              pow(r,5)*(-407.8425333707173 - 314.37502635182875*t + 0.5609391702828131*pow(t,2)) +
              r*(-0.1599506126960389 + 0.32576916425958746*t + 1.054371573642127*pow(t,2)) +
              pow(r,3)*(53.00757425014704 + 47.767873798856755*t + 8.477385842364002*pow(t,2)))/
              (pow(0.15247058823529414 + 0.8470588235294118*r - 1.7294117647058822*pow(r,2) +
              pow(r,3),2)*pow(pow(0.648 + 3.6*r - 7.35*pow(r,2) + 4.25*pow(r,3),2)*
              (-0.8 + 3.0*r - 1.25*pow(r,2) - t),0.3333333333333333)*(0.64 - 2.4*r + pow(r,2) + 0.8*t));
      dX_dt = -((2*(3 - 2.5*r)*(0.648 + 3.6*r - 7.35*pow(r,2) + 4.25*pow(r,3)) +
              (3.6 - 14.7*r + 12.75*pow(r,2))*(-0.8 + 3*r - 1.25*pow(r,2) - t))*((3 - 2.5*r)*pow(0.648 + 3.6*r - 7.35*pow(r,2) + 4.25*pow(r,3),2) +
              2*(3.6 - 14.7*r + 12.75*pow(r,2))*(0.648 + 3.6*r - 7.35*pow(r,2) + 4.25*pow(r,3))*(-0.8 + 3*r - 1.25*pow(r,2) - t)))/
              (3.0*pow(6,0.3333333333333333)*pow(pow(0.648 + 3.6*r - 7.35*pow(r,2) + 4.25*pow(r,3),2)*(-0.8 + 3*r - 1.25*pow(r,2) - t),1.3333333333333333)) +
              (3*(3 - 2.5*r)*(3.6 - 14.7*r + 12.75*pow(r,2)) - 5.*(0.648 + 3.6*r - 7.35*pow(r,2) + 4.25*pow(r,3)) +
              (-14.7 + 25.5*r)*(-0.8 + 3*r - 1.25*pow(r,2) - t))/(pow(6,0.3333333333333333)*pow(pow(0.648 + 3.6*r - 7.35*pow(r,2) +
              4.25*pow(r,3),2)*(-0.8 + 3*r - 1.25*pow(r,2) - t),0.3333333333333333));
      Y = 1.6509636244473134*pow((0.648 + 3.6*r - 7.35*pow(r,2) + 4.25*pow(r,3))*pow(-0.8 + 3*r - 1.25*pow(r,2) - t,2),0.3333333333333333);
      dY_dr = (-0.44377902225143756 - 143.16950180754043*pow(r,5) + 25.581337410056026*pow(r,6) + pow(r,3)*(-280.0859788874867 - 96.58137203016783*t) + 1.030201301655124*t +
              1.9811563493367759*pow(t,2) + pow(r,4)*(295.79764938014364 + 29.23581418292117*t) + r*(-16.351143736526197 - 34.93439029330516*t - 8.089721759791836*pow(t,2)) +
              pow(r,2)*(120.09109404229754 + 98.89272110439406*t + 7.016595403901082*pow(t,2)))/pow((0.648 + 3.6*r - 7.35*pow(r,2) + 4.25*pow(r,3))*pow(0.8 - 3*r + 1.25*pow(r,2) + t,2),0.6666666666666666);
      dY_dt = (5.847162836584234*(0.15247058823529414 + 0.8470588235294118*r - 1.7294117647058822*pow(r,2) + pow(r,3))*(0.64 - 2.4*r + pow(r,2) + 0.8*t))/
              pow((0.648 + 3.6*r - 7.35*pow(r,2) + 4.25*pow(r,3))*pow(0.8 - 3*r + 1.25*pow(r,2) + t,2),0.6666666666666666);
    }
    double H_Denom = X * Y * pow(pow(X,2) - pow(fp,2), 1.5);
    if(isnan(H_Denom))//If there is a negative number under the radical, throw an exception
    {
        printf("Error! Cannot Have a Negative Number Raised to the Power 3/2!\n");
        printf("r= %lf f'= %lf X= %lf\n", r, fp, X);
        return NULL;
    }
    else if(H_Denom == 0.0)//If the denominator = 0, throw an exception
    {
        printf("Error! Cannot Divide by 0!\n");
        printf("r= %lf f'= %lf X= %lf\n", r, fp, X);
        return NULL;
    }
    double H_Num = (2.0*pow(fp,3)*dY_dr) - (pow(X,3)*Y*dX_dt) + ((X*Y*fp)*(dX_dr + 2*fp*dX_dt)) - (2.0*pow(X,4)*dY_dt) -(pow(X,2)*((Y*fpp)+(2.0*fp)*(dY_dr - fp*dY_dt)));
    H = H_Num / H_Denom;
    printf("r: %.5lf fbar': %.15lf fbar'': %.15lf\n", r, fp - 1, fpp);
    double* H_ptr = &H;
    return H_ptr;//Return a pointer to the value of H
}

/*Function needed to expose the C function paraEvalFbar()

Python parameters: float variable x representing value of x

Python return: tuple representing (r, fbar(r) = t) evaluated at the input x value
value
*/
static PyObject* fbarParaEval(PyObject* self, PyObject* args)
{
  double x;
  if(!PyArg_ParseTuple(args, "d", &x))//If x is not a double, throw an exception
    return NULL;
  double* rfbar = paraEvalFbar(x);//Perform numerical evaluation
  PyObject* result = Py_BuildValue("dd", rfbar[0], rfbar[1]);//Build value to be returned to Python script
  free(rfbar);//Free memory alloctaed for the instantiation of the array
  return result;
}

/*Function needed to expose the C function explicitEvalFbar()

Python parameters: float variable r representing value of r

Python return: double that is the value of fbar at the input r value
*/
static PyObject* fbarExplicitEval(PyObject* self, PyObject* args)
{
  double r;//Value of parameter x(from Python)
  if(!PyArg_ParseTuple(args, "d", &r))//If x is not a double, throw an exception
    return NULL;
  double fbar = explicitEvalFbar(r);//Perform numerical evaluation
  return Py_BuildValue("d", fbar);//Build value to be returned to Python script
}

/*Function needed to expose the C function evalDerivativesPara()

Python paramters: float variable x representing value of x, float variable r representing value of r

Python return: tuple representing(fbar'(r), fbar''(r)) evaluated at the input x- and r- values
*/
static PyObject* calcDerivativesPara(PyObject* self, PyObject* args)
{
  double x, r;//Values of x and r, from Python
  if(!PyArg_ParseTuple(args, "dd", &x, &r))//If x is not a double, throw an exception
    return NULL;
  double* deriv = evalDerivativesPara(x, r);//Perform differentiation
  PyObject* derivPy = Py_BuildValue("dd", deriv[0], deriv[1]);//Build value to be returned to Python script
  free(deriv);//Free memory allocated during instantiation of array
  return derivPy;//Return 2-tuple to Python script
}

/*Function needed to expose the C function evalDerivativesExplicit()

Python paramters: float variable x representing value of x, float variable r representing value of r

Python return: tuple representing(fbar'(r), fbar''(r)) evaluated at the input x- and r- values
*/
static PyObject* calcDerivativesExplicit(PyObject* self, PyObject* args)
{
  double r;//Value of r, from Python
  if(!PyArg_ParseTuple(args, "d", &r))
    return NULL;//If r is not a double, throw an exception
  double* deriv = evalDerivativesExplicit(r);//Perform differentiation
  PyObject* derivPy = Py_BuildValue("dd", deriv[0], deriv[1]);//Build value to be returned to Python script
  free(deriv);//Free memory allocated during instantiation of array
  return derivPy;//Return 2-tuple to Python script
}

/*Function needed to expose the C function calculateHS()

Python paramters: float variables representing values of r, fbar, fbar'(r), and fbar''(r) = f''(r)

Python return: float representing value of H caluclated using the paramter values
*/
static PyObject* calcHS1(PyObject* self, PyObject* args)
{
  double r, fbar, fbarp, fpp;//Value of r, f(r) = t, f'(r), and f''(r) from Python
  if(!PyArg_ParseTuple(args, "dddd", &r, &fbar, &fbarp, &fpp))//If any of the inputs are not doubles, throw an exception
    return NULL;
  double* H_ptr = calculateHS(r, fbar, fbarp, fpp);//Calculate H
  return Py_BuildValue("d", *H_ptr);//Return what H_ptr points to
}

/*Function needed to expose the C function calculateHGF()

Python paramters: float variables representing values of r, t, f'(r), and f''(r)

Python return: float representing value of H caluclated using the paramter values
*/
static PyObject* calcHGF1(PyObject* self, PyObject* args)
{
  double r, t, fp, fpp;//Value of r, f(r) = t, f'(r), and f''(r) from Python
  if(!PyArg_ParseTuple(args, "dddd", &r, &t, &fp, &fpp))//If any of the inputs are not doubles, throw an exception
    return NULL;
  double* H_ptr = calculateHGF(r, t, fp, fpp);//Calculate H
  return Py_BuildValue("d", *H_ptr);//Return what H_ptr points to
}

//Mapping of the names of Python methods to C functions
static PyMethodDef fulmineMethods[] = {
      {"peFbar", fbarParaEval, METH_VARARGS, evalFbarPara_docstring},//Mapping to paraEvalFbar()
      {"eeFbar", fbarExplicitEval, METH_VARARGS, evalFbarExplicit_docstring},//Mapping to explicitEvalFbar()
      {"peCalcDeriv", calcDerivativesPara, METH_VARARGS, derivPara_docstring},//Mapping to evalDerivativesPara()
      {"eeCalcDeriv", calcDerivativesExplicit, METH_VARARGS, derivExplicit_docstring},//Mapping to evalDerivativesExplicit()
      {"calcHS", calcHS1, METH_VARARGS, calcHS_docstring},//Mapping to calculateHS()
      {"calcHGF", calcHGF1, METH_VARARGS, calcHGF_docstring},//Mapping to calculateHGF()
      {NULL, NULL, 0, NULL}//Needed for formatting
};

//Module definition structure
static struct PyModuleDef fulmineModule = {
    PyModuleDef_HEAD_INIT,
    "fulmine",             //Name of module
    module_docstring,      //Documentation
    -1,
    fulmineMethods
};

//Function to initialize the module
PyMODINIT_FUNC PyInit_fulmine(void)
{
    return PyModule_Create(&fulmineModule);
}

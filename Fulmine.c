/*Eric Gelphman
  University of California San Diego(UCSD)
  Department of Mathematics
  Jacobs School of Engineering Department of Electrical and Computer Engineering(ECE)
  September 17, 2017

  Fulmine, a C module to perform arithmetic operations for fGenerator and KrazyGlue.
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
static char calcH_docstring[] = "Function to calculate the mean curvature H";
static char calcHLNr_docstring[] = "Function to calculate the mean curvature H when r is a large negative number";

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
    //Evaluation of r in the parametric region
    if(x < 1.0e-16)//10th order Taylor series approximation for r for very small x
    {
      r = 5.333333333333333 - (5*x)/3.0 + (19*pow(x,2))/12.0 + (17*pow(x,3))/12.0 + (131*pow(x,4))/96.0 + (323*pow(x,5))/240.0 + (257*pow(x,6))/192.0 +
          (1795*pow(x,7))/1344.0 + (4099*pow(x,8))/3072.0 + (3073*pow(x,9))/2304.0 + (20483*pow(x,10))/15360.0 - 2.0*log(2.0) + 2.0*log(x);
      //Note that log(x) = ln(x)
    }
    else//Regular evluation for r
      r = 4.0/(3.0*a) + 4.0*a - 4.0*atanh(a);
    result[0] = r;
    if(r < -60.0)//10th order Taylor series approximation in parametric region for large negative r
    {
      fbar = -2.0*x - (5*pow(x,2))/2.0 - (35*pow(x,3))/12.0 - (105*pow(x,4))/32.0 - (231*pow(x,5))/64.0 - (1001*pow(x,6))/256.0 -
             (2145*pow(x,7))/512.0 - (36465*pow(x,8))/8192.0 - (230945*pow(x,9))/49152.0 - (323323*pow(x,10))/65536.0;
    }
    else//Regular evaluation for fbar
      fbar = 4.0/3.0 - 4.0/(3.0*pow(a,1.5));
    result[1] = fbar;
    return result;
}

/*Function to evaluate fbar in the glue(-1.0 < r < -0.6) and linear(r >= -0.6) regions,
where it is explicitly defined in terms of the independent variable r

Parameter r is the value of r that fbar is to be evaluated at

Returns a double which is the value of fbar for the input value of r
*/
double explicitEvalFbar(double r)
{
    double result = 0.0;
    if(r > -1.0 && r < -0.6)//Evaluation in glue region
      result = 5.4215319*pow(r,3) + 11.9073083*pow(r,2) + 7.4335155*r + 0.7445292;
    else if(r >= -0.6)//Evaluation in linear region
      result = -1.0*r - 1.2;
    return result;
}

/*Function to calculate the first and second order derivatives of fbar in
the parametric region

For r <= -100.0, taylor series of derivative is used for the numerical evaluation,
otherwise the derivative is evaluated using its explicit formula

Parameter x is the value of x the formula is to be evaluated at, r is the corresponding r-value

Returns a pointer to an array which stores the values of the first and second derivatives
*/
double* evalDerivativesPara(double x, double r)
{
    double* deriv = malloc(sizeof(double)*2);
    double d1, d2;//fbar', fbar''
    if(r <= -60.0)//10th-order Taylor series expansion
    {
      d1 = -x - 3.0*pow(x,2) - (25*pow(x,3))/8.0 + (37*pow(x,4))/8.0 + (2817*pow(x,5))/128.0 + (1891*pow(x,6))/64.0 - (22129*pow(x,7))/1024.0 -
      (165207*pow(x,8))/1024.0 - (8666971*pow(x,9))/32768.0 + (225585*pow(x,10))/4096.0;
      d2 = -x/2.0 - (13*pow(x,2))/4.0 - (79*pow(x,3))/16.0 + (575*pow(x,4))/32.0 + (24453*pow(x,5))/256.0 + (69255*pow(x,6))/512.0 -
           (580671*pow(x,7))/2048.0 - (6626277*pow(x,8))/4096.0 - (164430675*pow(x,9))/65536.0 + (381990279*pow(x,10))/131072.0;
    }
    else//Explicit evaluation
    {
      d1 = -2/((-4 + 2/pow(1 - x,2.5) - 4/((-2 + x)*x))*pow(1 - x,2.5));
      d2 = -((-2 + x)*pow(-1 + x,6)*x*(-4 - 10*x + 5*pow(x,2)))/(2.0*pow(2*sqrt(1 - x) + (2 - 8*sqrt(1 - x))*x +
           (-1 + 12*sqrt(1 - x))*pow(x,2) - 8*sqrt(1 - x)*pow(x,3) + 2*sqrt(1 - x)*pow(x,4),3));
    }
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
    if(r > -1.0 && r < -0.6)//Glue region
    {
      dfbar_1 = 16.2645957*pow(r,2) + 23.8146166*r + 7.4335155;
      dfbar_2 = 32.5291914*r + 23.8146166;
    }
    else//Linear region (r >= -0.6)
    {
      dfbar_1 = -1.0;
      dfbar_2 = 0.0;
    }
    deriv[0] = dfbar_1;
    deriv[1] = dfbar_2;
    return deriv;
}

/*Function to calculate H for when r is a large negative number

Paramters: r, t, fp, fpp are the values of r, t, f'(r), and f''(r) from fGenerator(fbar converted to f)

Return: Pointer to a double that stores the value of the mean curvature H
*/
double* calculateHLNr(double r, double fbar, double fbarp, double fpp)
{
  double X, Xbar, Y, a, b, z, H_Num, H_Denom;
  static double H;
  z = -fbar + 4.0/3.0;//z = r - t, t = fbar + r - 4/3
  a = pow(z, 1.0/3.0);
  b = pow(4.0/3.0, 1.0/3.0);
  //Express Xbar = X - 1 as a 10th-order Taylor series centered around (r - t)^(1/3) = (4/3)^(1/3)
  Xbar = -((pow(3,0.3333333333333333)*(a-b))/pow(2,0.6666666666666666)) + (pow(3,0.6666666666666666)*pow(a-b,2))/(2.0*pow(2,0.3333333333333333)) -
         (3*pow(a-b,3))/4.0 + (3*pow(3,0.3333333333333333)*pow(a-b,4))/(4.0*pow(2,0.6666666666666666)) -
         (3*pow(3,0.6666666666666666)*pow(a-b,5))/(8.0*pow(2,0.3333333333333333)) + (9*pow(a-b,6))/16.0 -
         (9*pow(3,0.3333333333333333)*pow(a-b,7))/(16.0*pow(2,0.6666666666666666)) + (9*pow(3,0.6666666666666666)*pow(a-b,8))/(32.0*pow(2,0.3333333333333333)) -
         (27*pow(a-b,9))/64.0 + (27*pow(3,0.3333333333333333)*pow(a-b,10))/(64.0*pow(2,0.6666666666666666));
  Y = pow(4.50, 1.0/3.0)*pow(z, 2.0/3.0);
  X = pow(2,0.6666666666666666)/(pow(3,0.3333333333333333)*pow(z,0.3333333333333333));
  H_Denom = X * Y * pow(pow(X,2) - pow(fbar,2), 1.5);
  if(H_Denom == 0.0)//If the denominator = 0, throw an exception
  {
    printf("Error! Cannot Divide by 0!\n");
    return NULL;
  }
  //double aleph = 2.0 + (-fbar)/2.0 + pow(-fbar,2)/4.0 - (7*pow(-fbar,3))/48.0 + (35*pow(-fbar,4))/384.0 - (91*pow(-fbar,5))/1536.0;
  double alpha = 11.0 - (9*(-fbar))/4.0 + (21*pow(-fbar,2))/16.0 - (55*pow(-fbar,3))/64.0 + (151*pow(-fbar,4))/256.0 - (425*pow(-fbar,5))/1024.0;
  double beta = 19.5 - (21*(-fbar))/8.0 + (51*pow(-fbar,2))/32.0 - (137*pow(-fbar,3))/128.0 + (383*pow(-fbar,4))/512.0 -
                (1093*pow(-fbar,5))/2048.0;
  double gamma = -0.5 + (5*(-fbar))/8.0 - (5*pow(-fbar,2))/8.0 + (55*pow(-fbar,3))/96.0 - (385*pow(-fbar,4))/768.0 + (1309*pow(-fbar,5))/3072.0 -
                 (6545*pow(-fbar,6))/18432.0;
  double delta = 2.0 - (5*(-fbar))/2.0 + (5*pow(-fbar,2))/2.0 - (55*pow(-fbar,3))/24.0 + (385*pow(-fbar,4))/192.0 - (1309*pow(-fbar,5))/768.0 +
                 (6545*pow(-fbar,6))/4608.0;
  double epsilon = 0.5 - (3*(-fbar))/8.0 + (9*pow(-fbar,2))/32.0 - (27*pow(-fbar,3))/128.0 + (81*pow(-fbar,4))/512.0 - (243*pow(-fbar,5))/2048.0 +
                   (729*pow(-fbar,6))/8192.0;
  double zeta = 2.0 + fbar/2.0 + pow(-fbar,2)/4.0 - (7*pow(-fbar,3))/48.0 + (35*pow(-fbar,4))/384.0 - (91*pow(-fbar,5))/1536.0 +
                (91*pow(-fbar,6))/2304.0;
  double eta = -4.0 + 3*(-fbar) - (9*pow(-fbar,2))/4.0 + (27*pow(-fbar,3))/16.0 - (81*pow(-fbar,4))/64.0 + (243*pow(-fbar,5))/256.0 -
               (729*pow(-fbar,6))/1024.0;
  H_Num = 2.0*X*pow(fbarp, 3) + alpha*pow(fbarp, 2) + beta*fbarp + gamma + delta + epsilon + zeta + eta - 2.0*fpp;
  H = H_Num / H_Denom;
  double* H_Ptr = &H;
  return H_Ptr;
}

/* Function to calculate the numerical values of X, Y, and their first order partial
derivatives

Paramters: r, t, fp, fpp are the values of r, t, f'(r), and f''(r) from fGenerator(fbar converted to f)

Return: Pointer to a double that stores the value of the mean curvature H
*/
double* calculateH(double r, double t, double fp, double fpp)
{
    double X, Y, dX_dr, dX_dt, dY_dr, dY_dt;
    static double H;
    if(r <= 0.8)//Schwarzchild Spacetime M(r) = 1, t0(r) = r, applies for all r <= 0.8
    {
      X = pow(2,0.6666666666666666)/(pow(3,0.3333333333333333)*pow(r - t,0.3333333333333333));
      //printf("X: %lf\n", X);
      dX_dr = -pow(2,0.6666666666666666)/(3.0*pow(3,0.3333333333333333)*pow(r - t,1.3333333333333333));
      dX_dt = -dX_dr;
      Y = 1.6509636244473134*pow(r - t, 2.0/3.0);
      dY_dr = X;
      dY_dt = -X;//dY/dt = -dY/dr = -X
    }
    else if(r >= 1.2)//Friedman Spacetime M(r) = r^3, t0(r) = 1, applies for all r >= 1.2
    {
      X = (pow(3,0.6666666666666666)*pow(r,2)*(1 - t))/(pow(2,0.3333333333333333)*pow(pow(r,6)*(1.0 - t),0.3333333333333333));
      dX_dr = 0.0;
      dX_dt = -((pow(2,0.6666666666666666)*pow(r,2))/(pow(3,0.3333333333333333)*pow(-(pow(r,6)*(-1.0 + t)),0.3333333333333333)));
      Y = 1.6509636244473134*r*pow(1 - t, 2.0/3);
      dY_dr = 1.6509636244473134*pow(1 - t, 2.0/3);
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
    if(fabs(X) < fabs(fp))//If there is a negative number under the radical, throw an exception
    {
        printf("Error! Cannot Have a Negative Number Raised to the Power 3/2!\n");
        printf("r: %lf f': %lf X: %lf\n", r, fp, X);
        return NULL;
    }
    double H_Denom = X * Y * pow(pow(X,2) - pow(fp,2), 1.5);
    if(H_Denom == 0.0)//If the denominator = 0, throw an exception
    {
        printf("Error! Cannot Divide by 0!\n");
        printf("r: %lf f': %lf X: %lf\n", r, fp, X);
        return NULL;
    }
    double H_Num = (2.0*pow(fp,3)*dY_dr) - (pow(X,3)*Y*dX_dt) + ((X*Y*fp)*(dX_dr + 2*fp*dX_dt)) - (2*pow(X,4)*dY_dt) -
                   (pow(X,2)*((Y*fpp)+(2*fp)*(dY_dr - fp*dY_dt)));
    H = H_Num / H_Denom;
    //printf("r: %lf H: %lf\n", r, H);
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

/*Function needed to expose the C function calculateHLNr()

Python paramters: float variables representing values of r, fbar, fbar'(r), and fbar''(r) = f''(r)

Python return: float representing value of H caluclated using the paramter values
*/
static PyObject* calcHLNr1(PyObject* self, PyObject* args)
{
  double r, fbar, fbarp, fpp;//Value of r, f(r) = t, f'(r), and f''(r) from Python
  if(!PyArg_ParseTuple(args, "dddd", &r, &fbar, &fbarp, &fpp))//If any of the inputs are not doubles, throw an exception
    return NULL;
  double* H_ptr = calculateHLNr(r, fbar, fbarp, fpp);//Calculate H
  return Py_BuildValue("d", *H_ptr);//Return what H_ptr points to
}

/*Function needed to expose the C function calculateH()

Python paramters: float variables representing values of r, t, f'(r), and f''(r)

Python return: float representing value of H caluclated using the paramter values
*/
static PyObject* calcH1(PyObject* self, PyObject* args)
{
  double r, t, fp, fpp;//Value of r, f(r) = t, f'(r), and f''(r) from Python
  if(!PyArg_ParseTuple(args, "dddd", &r, &t, &fp, &fpp))//If any of the inputs are not doubles, throw an exception
    return NULL;
  double* H_ptr = calculateH(r, t, fp, fpp);//Calculate H
  return Py_BuildValue("d", *H_ptr);//Return what H_ptr points to
}

//Mapping of the names of Python methods to C functions
static PyMethodDef fulmineMethods[] = {
      {"peFbar", fbarParaEval, METH_VARARGS, evalFbarPara_docstring},//Mapping to paraEvalFbar()
      {"eeFbar", fbarExplicitEval, METH_VARARGS, evalFbarExplicit_docstring},//Mapping to explicitEvalFbar()
      {"peCalcDeriv", calcDerivativesPara, METH_VARARGS, derivPara_docstring},//Mapping to evalDerivativesPara()
      {"eeCalcDeriv", calcDerivativesExplicit, METH_VARARGS, derivExplicit_docstring},//Mapping to evalDerivativesExplicit()
      {"calcHLNr", calcHLNr1, METH_VARARGS, calcHLNr_docstring},//Mapping to calculateHLNr()
      {"calcH", calcH1, METH_VARARGS, calcH_docstring},//Mapping to calculateH()
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

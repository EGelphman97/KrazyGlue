/*Eric Gelphman
  University of California San Diego(UCSD)
  Department of Mathematics
  Jacobs School of Engineering Department of Electrical and Computer Engineering(ECE)
  September 10, 2017

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
    if(x < 5.551116e-17)//Linear approximation for r for very small x
      r = 2*log(x) + 16/3 - log(4) - x;//Note that log(x) = ln(x)
    else//Regular evluation for r
      r = 4.0/(3*a) + 4*a - 4*atanh(a);
    result[0] = r;
    fbar = ((double)4.0)/3.0 - 4.0/(3.0*pow(a,1.5));
    /*
    if(r < -60.0)//10th Order Taylor Series Approximation in parametric region for large negative r
      fbar = -2.0*x + 2.5*pow(x,2) - 2.91667*pow(x,3) - 3.28125*pow(x,4) - 3.60938*pow(x,5) - 3.91016*pow(x,6) - 4.18945*pow(x,7) - 4.45129*pow(x,8) - 4.70117*pow(x,9) - 4.93352*pow(x,10);
      */
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
      result = -13.49319*pow(r,3) - 33.45508*pow(r,2) - 26.57345*r - 6.81477;
    else if(r >= -0.6)//Evaluation in linear region
      result = -r - 0.6;
    return result;
}

/*Function to calculate the first and second order derivatives of fbar in
the parametric region

Though it is possible to calculate the derivatives explicitly, doing so would likely
induce large numerical errors due to how the atanh() function is defined.
Thus, derivatives are calculated by using 10th-order series expansions of both r and fbar
and then applying the parametric differentiation formulas. This results in an
approximation of a very high degree.

Parameter x is the value of x the formula is to be evaluated at

Returns a pointer to an array which stores the values of the first and second derivatives
*/
double* evalDerivativesPara(double x)
{
    double* deriv = malloc(sizeof(double)*2);
    //fbar'
    deriv[0] = (-3.0*x*(65536 + 163840*x + 286720*pow(x,2) + 430080*pow(x,3) + 591360*pow(x,4) +
      768768*pow(x,5) + 960960*pow(x,6) + 1166880*pow(x,7) + 1385670*pow(x,8) +
      1616615*pow(x,9))) / (64.0*(3072 - 2560*x + 4864*pow(x,2) + 6528*pow(x,3) + 8384*pow(x,4) +
      10336*pow(x,5) + 12336*pow(x,6) + 14360*pow(x,7) + 16396*pow(x,8) +
      18438*pow(x,9) + 20483*pow(x,10)));
    //fbar''
    deriv[1] = (-576*x*(25165824 + 125829120*x + 238026752*pow(x,2) + 370147328*pow(x,3) +
      557121536*pow(x,4) + 855179264*pow(x,5) + 1341808640*pow(x,6) +
      2117083136*pow(x,7) + 3305365760*pow(x,8) + 5057066240*pow(x,9) -
      302526208*pow(x,10) + 9052516800*pow(x,11) + 7388817920*pow(x,12) +
      5715151248*pow(x,13) + 4092533456*pow(x,14) + 2642954028*pow(x,15) +
      1474781880*pow(x,16) + 651204125*pow(x,17) + 178058595*pow(x,18))) /
      pow(3072 - 2560*x + 4864*pow(x,2) + 6528*pow(x,3) + 8384*pow(x,4) +
      10336*pow(x,5) + 12336*pow(x,6) + 14360*pow(x,7) + 16396*pow(x,8) +
      18438*pow(x,9) + 20483*pow(x,10),3);
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
    double dfbar_1 = 0.0;
    double dfbar_2 = 0.0;//First, and second derivatives, respectively
    double* deriv = malloc(sizeof(double)*2);
    if(r > -1.0 && r < -0.6)//Glue region
    {
      dfbar_1 = -40.47957*pow(r,2) - 66.91016*r -26.57345;
      dfbar_2 = -80.95914*r - 66.91016;
    }
    else if(r >= -0.6)//Linear region
    {
      dfbar_1 = -1.0;
      dfbar_2 = 0.0;
    }
    deriv[0] = dfbar_1;
    deriv[1] = dfbar_2;
    return deriv;
}

/* Function to calculate the numerical values of X, Y, and their first order partial
derivatives

Paramters: r, t, fp, fpp are the values of r, t, f'(r), and f''(r) from fGenerator(fbar converted to f)

Return: Pointer to a double that stores the value of the mean curvature H
*/
double* calculateH(double r, double t, double fp, double fpp)
{
    double X, Y, dX_dr, dX_dt, dY_dr, dY_dt, H, H_Denom;
    if(r <= 0.8)//Schwarzchild Region
    {
      X = pow(2,0.6666666666666666)/(pow(3,0.3333333333333333)*pow(r - t,0.3333333333333333));
      dX_dr = -pow(2,0.6666666666666666)/(3.0*pow(3,0.3333333333333333)*pow(r - t,1.3333333333333333));
      dX_dt = -dX_dr;
      Y = 1.6509636244473134*pow(pow(r - t,2),0.3333333333333333);
      dY_dr = 1.10064/pow(r - t,0.3333333333333333);
      dY_dt = -dY_dr;
    }
    else if(r > 0.8 && r < 1.2)//Glue Region
    {
      X = (pow(3,0.6666666666666666)*pow(r,2)*(1 - t))/(pow(2,0.3333333333333333)*pow(pow(r,6)*(1 - t),0.3333333333333333));
      dX_dr = 0.0;
      dX_dt = -((pow(2,0.6666666666666666)*pow(r,2))/(pow(3,0.3333333333333333)*pow(-(pow(r,6)*(-1 + t)),0.3333333333333333)));
      Y = pow(2,0.6666666666666666)/(pow(3,0.3333333333333333)*pow(r - t,0.3333333333333333));
      dY_dr = (1.6509636244473134*pow(r,2)*pow(1 - t,2))/pow(pow(r,3)*pow(1 - t,2),0.6666666666666666);
      dY_dt = (1.6509636244473134*pow(r,2)*pow(1 - t,2))/pow(pow(r,3)*pow(1 - t,2),0.6666666666666666);
    }
    else//r >= 1.2 Friedman Region
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
    if(fabs(X) < fabs(fp))
        return NULL;
    H_Denom = X * Y * pow(pow(X,2) - pow(fp,2), 1.5);
    if(H_Denom == 0.0)
        return NULL;
    H = ((2.0*pow(fp,3)*dY_dr) - (pow(X,3)*Y*dX_dt) + ((X*Y*fp)*(dX_dr + 2*fp*dX_dt)) - (2*pow(X,4)*dY_dt) - (pow(X,2)*((Y*fpp)+(2*fp)*(dY_dr - fp*dY_dt)))) / H_Denom;
    double* H_ptr = &H;
    return H_ptr;
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
  double x;//Values of x and r, from Python
  if(!PyArg_ParseTuple(args, "d", &x))//If x is not a double, throw an exception
    return NULL;
  double* deriv = evalDerivativesPara(x);//Perform differentiation
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

/*Function needed to expose the C function calculateH()

Python paramters: float variables representing values of r, t, f'(r), and f''(r)

Python return: float representing value of H caluclated using the paramter values
*/
static PyObject* calcH1(PyObject* self, PyObject* args)
{
  double r, t, fp, fpp;//Value of r, f(r) = t, f'(r), and f''(r) from Python
  if(!PyArg_ParseTuple(args, "dddd", &r, &t, &fp, &fpp))//If r and t are not doubles, throw an exception
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

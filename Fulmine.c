/*Eric Gelphman
  University of California San Diego(UCSD)
  Department of Mathematics
  Jacobs School of Engineering Department of Electrical and Computer Engineering(ECE)
  September 7, 2017

  Fulmine, a C subroutine to perform numerical differentiation and evaluation
  for fGenerator and KrazyGlue.
*/

#include <Python.h>//Library needed to interface with Python
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

/*Function to evaluate r and fbar in the parametric region(r <= -1.0), where both
r and fbar(r) = t are defined in terms of the parameter x

Parameter x is the value of x that r and fbar are to evaluated at

Returns a pointer to an array which stores the values of r and fbar
*/
double* paraEvalFbar(double x)
{

    double a,r,fbar;
    double result[2];
    double* resultPtr;
    resultPtr = result;
    //Evaluation of r in the parametric region
    if(x < 5.551116e-17)//Linear approximation for r for very small x
    { r = 2*log(x) + 16/3 - log(4) - x; }
    else//Regular evluation for r
    { r = 4.0/(3*a) + 4*a - 4*atanh(a); }
    result[0] = r;
    if(r < -60.0)//10th Order Taylor Series Approximation in parametric region for large negative r
    { fbar = -2*x + 2.5*pow(x,2) - 2.91667*pow(x,3) - 3.28125*pow(x,4) - 3.60938*pow(x,5) - 3.91016*pow(x,6) - 4.18945*pow(x,7) - 4.45129*pow(x,8) - 4.70117*pow(x,9) - 4.93352*pow(x,10); }
    else if(r >= -60.0 && r <= -1.0)//Regular evaluation in parametric region
    { fbar = 4.0/3 - 4.0/(3*pow(a,1.5)); }
    else if(r > -1.0 && r < -0.6)//Evaluation in glue region
    { fbar = -13.49319*pow(r,3) - 33.45508*pow(r,2) - 26.57345*r - 6.81477; }
    else// r >= -0.6, in linear region
    { fbar = -r - 0.6; }
    result[1] = fbar;
    return resultPtr;
}

/*Function to evaluate fbar in the glue(-1.0 < r < -0.6) and linear(r >= -0.6) regions,
where it is explicitly defined in terms of the independent variable r

Parameter r is the value of r that fbar is to be evaluated at

Returns a double which is the value of fbar for the input value of r
*/
double explicitEvalFbar(double r)
{
    if(r > -1.0 && r < -0.6)//Evaluation in glue region
    { return -13.49319*pow(r,3) - 33.45508*pow(r,2) - 26.57345*r - 6.81477; }
    else if(r >= -0.6)//Evaluation in linear region
    { return -r - 0.6; }
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
    double deriv[2];
    double* derivPtr;
    derivPtr = deriv;
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
    return derivPtr;
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
    double deriv[2];
    double* derivPtr;
    derivPtr = deriv;
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
    return derivPtr;
}

/*Function needed to expose the C function paraEvalFbar()

Python parameters: float variable x representing value of x

Python return: tuple representing (r, fbar(r) = t) evaluated at the input x value
value
*/
static PyObject* fbarParaEval(PyObject* self, PyObject* args)
{
  double x;//Value of parameter x(from Python)
  if(!PyArgParseTuple(args, "d", &x))//If x is not a double, throw an exception
  { return NULL; }
  double* rfbar = paraEvalFbar(x);//Perform numerical evaluation
  return PyBuildValue("dd", *rfbar, *(rfbar + 1));//Build value to be returned to Python script
}

/*Function needed to expose the C function explicitEvalFbar()

Python parameters: float variable r representing value of r

Python return: double that is the value of fbar at the input r value
*/
static PyObject* fbarExplicitEval(PyObject* self, PyObject* args)
{
  double r;//Value of parameter x(from Python)
  if(!PyArgParseTuple(args, "d", &r))//If x is not a double, throw an exception
  { return NULL; }
  double fbar = explicitEvalFbar(r);//Perform numerical evaluation
  return PyBuildValue("d", fbar);//Build value to be returned to Python script
}

/*Function needed to expose the C function evalDerivativesPara()

Python paramters: float variable x representing value of x, float variable r representing value of r

Python return: tuple representing(fbar'(r), fbar''(r)) evaluated at the input x- and r- values
*/
static PyObject* calcDerivativesPara(PyObject* self, PyObject* args)
{
  double x;//Values of x and r, from Python
  if(!PyArgParseTuple(args, "d", &x))//If x is not a double, throw an exception
  { return NULL; }
  double* deriv = evalDerivativesPara(x);//Perform differentiation
  return PyBuildValue("dd", *deriv, *(deriv + 1));//Build value to be returned to Python script
}

/*Function needed to expose the C function evalDerivativesExplicit()

Python paramters: float variable x representing value of x, float variable r representing value of r

Python return: tuple representing(fbar'(r), fbar''(r)) evaluated at the input x- and r- values
*/
static PyObject* calcDerivativesExplicit(PyObject* self, PyObject* args)
{
  double r;//Values of x and r, from Python
  if(!PyArgParseTuple(args, "d", &r))//If r is not a double, throw an exception
  { return NULL; }
  double* deriv = evalDerivativesExplicit(r);//Perform differentiation
  return PyBuildValue("dd", *deriv, *(deriv + 1));//Build value to be returned to Python script
}

//Mapping of the names of Python methods to C functions
static PyMethodDef moduleMethods[] = {
      {"peFbar", fbarParaEval, METH_VARARGS, evalFbarPara_docstring},//Mapping to paraEvalFbar()
      {"eeFbar", fbarExplicitEval, METH_VARARGS, evalFbarExplicit_docstring},//Mapping to explicitEvalFbar()
      {"peCalcDeriv", calcDerivativesPara, METH_VARARGS, derivPara_docstring},//Mapping to evalDerivativesPara()
      {"eeCalcDeriv", calcDerivativesExplicit, METH_VARARGS, derivExplicit_docstring},//Mapping to evalDerivativesExplicit()
      {NULL, NULL, 0, NULL}//Needed for formatting
};

//Function to initialize the module
PyMODINIT_FUNC init_Fulmine(void)
{
    PyObject* fulmine = Py_InitModule3("Fulmine", moduleMethods, module_docstring);
    if(fulmine = NULL)
    { return; }
}

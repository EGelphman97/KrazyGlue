#!/usr/bin/env wolframscript
(* ::Package:: *)

ClearAll[format1]
format1 = RawBoxes[Replace[ToBoxes@#, InterpretationBox[a_, b_, c___] :>
      With[{aa = StringReplace[a, 
            { "Sqrt" -> "sqrt", "Power(E," -> "exp(", "Power" -> "pow", "ProductLog" -> "LambertW", "(-Power(w,2) + Power(z,2))/E" -> "a", "E" -> "exp(1.0)"}]}, aa], 
      {0, Infinity}]] &;
U = w-z;
V = w+z;
t = 2*Log[-V/U];
r = 2*ProductLog[-(U*V)/E]
abs = (Sqrt[2]+Sqrt[r])/(Sqrt[2]-Sqrt[r]);
mabs = -abs;
tau = t + 2*Sqrt[2*r] - 2*Log[abs];
mtau = t + 2*Sqrt[w*r] - 2*Log[mabs];
rho = tau + Sqrt[(2/9)*r^3];
mrho = mtau + Sqrt[(2/9)*r^3];
drhodz = D[rho,z];
dmrhodz = D[mrho,z];
drhodw = D[rho,w];
dmrhodw = D[mrho,w];
dtaudz = D[tau,z];
dmtaudz = D[mtau,z];
dtaudw = D[tau,w];
dmtaudw = D[mtau,w];
Print[format1@CForm[Simplify[drhodz]]];
Print[format1@CForm[Simplify[dmrhodz]]];
Print[format1@CForm[Simplify[drhodw]]];
Print[format1@CForm[Simplify[dmrhodw]]];
Print[format1@CForm[Simplify[dtaudz]]];
Print[format1@CForm[Simplify[dmtaudz]]];
Print[format1@CForm[Simplify[dtaudw]]];
Print[format1@CForm[Simplify[dmtaudw]]];








      
      




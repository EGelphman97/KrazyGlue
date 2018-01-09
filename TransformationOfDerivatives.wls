#!/usr/bin/env wolframscript
(* ::Package:: *)

ClearAll[format1]
format1 = RawBoxes[Replace[ToBoxes@#, InterpretationBox[a_, b_, c___] :>
      With[{aa = StringReplace[a, 
            { "Sqrt" -> "sqrt", "Power(E," -> "exp(", "Power" -> "pow"}]}, aa], 
      {0, Infinity}]] &;
r = (4.5^(1/3))*(rho - tau)^(2/3);
a = r + 2*Log[Abs[(r/2)-1]];
et4 = Exp[(tau/4)-Sqrt[r/2]]*Sqrt[Abs[(Sqrt[2]+Sqrt[r])/(Sqrt[2]-Sqrt[r])]];
et4n = Exp[(tau/4)-Sqrt[r/2]]*Sqrt[-((Sqrt[2]+Sqrt[r])/(Sqrt[2]-Sqrt[r]))];
et4p = Exp[(tau/4)-Sqrt[r/2]]*Sqrt[(Sqrt[2]+Sqrt[r])/(Sqrt[2]-Sqrt[r])];
U = -Exp[a/4]/et4;
V = Exp[a/4]*et4;
z = (V-U)/2;
w = (V+U)/2;
Print[z/.{rho->0.1,tau->-0.408786}];
Print[w/.{rho->0.7,tau->-1.65}];
dwdrho = D[w, rho];
dzdrho = D[z,rho];
dwdz = dwdrho/dzdrho;
Print[dwdz/.{rho->0.7,tau->-1.65}];
d2wdz2 = D[dwdz,rho]/dzdrho;
Print[d2wdz2/.{rho->0.7,tau->-1.65}];






 

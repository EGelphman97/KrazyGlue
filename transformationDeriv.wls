#!/usr/bin/env wolframscript
(* ::Package:: *)

ClearAll[format1]
format1 = RawBoxes[Replace[ToBoxes@#, InterpretationBox[a_, b_, c___] :>
      With[{aa = StringReplace[a, 
            { "Sqrt" -> "sqrt", "Power(E," -> "exp(", "Power" -> "pow"}]}, aa], 
      {0, Infinity}]] &;
r = (4.5^(1/3))*(rho - tau)^(2/3);
a = r + 2*Log[(r/2)-1];
et4 = Exp[(tau/4)-Sqrt[r/2]]*Sqrt[Abs[(Sqrt[2]+Sqrt[r])/(Sqrt[2]-Sqrt[r])]];
U = -Exp[a/4]/et4;
V = Exp[a/4]*et4;
z = (V-U)/2;
w = (V+U)/2;
dzdp = D[z, rho];
dz2dp2 = D[dz_dp, rho];
dz2dpdT = D[dz_dp, tau];
dwdp = D[w, rho];
dw2dp2 = D[dw_dp, rho];
dw2dpdT = D[dw_dp, tau];
dwdT = D[w, tau];
dw2dT2 = D[dw_dT, tau];
dzdT = D[z, tau];
dz2dT2 = D[dz_dT, tau];
dgdp = (d1f*dzdp-dwdp)/(dwdT-d1f*dzdT);
dgdp = Simplify[dgdp];
Print["****************************dg/dp:******************************************"];
Print[format1@CForm[dgdp]];
d2gdp2 = (d2f*(dzdp+dzdT*dgdp)+d1f*(dz2dp2+(dz2dpdT*dgdp)+dgdp*(dz2dpdT+dz2dT2*dgdp))-dw2dp2+(dw2dpdT*dgdp)+(dw2dpdT+dw2dT2*dgdp)*dgdp)/(dwdT-(d1f*dzdT));
d2gdp2 = Simplify[d2gdp2];
Print["****************************d^2g/dp^2:******************************************"];
Print[format1@CForm[d2gdp2]];




 

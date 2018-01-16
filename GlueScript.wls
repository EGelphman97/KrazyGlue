#!/usr/bin/env wolframscript
(* ::Package:: *)

ClearAll[format1]
format1 = RawBoxes[Replace[ToBoxes@#, InterpretationBox[a_, b_, c___] :>
      With[{aa = StringReplace[a, 
            { "Sqrt" -> "sqrt", "Power(E," -> "exp(", "Power" -> "pow"}]}, aa], 
      {0, Infinity}]] &;
z1=0.4;
w1=0.0;
fp1=0.0;
fpp1=0.0;
z2=1.691486380082176;
w2=-0.7866899404272474;
fp2=-0.9288342771929007;
fpp2=-0.01601146279759889;
b={w1,fp1,fpp1,w2,fp2,fpp2};
A={{z1^5,z1^4,z1^3,z1^2,z1,1},{5*z1^4,4*z1^3,3*z1^2,2*z1,1,0},{20*z1^3,12*z1^2,6*z1,2,0,0},{z2^5,z2^4,z2^3,z2^2,z2,1},{5*z2^4,4*z2^3,3*z2^2,2*z2,1,0},{20*z2^3,12*z2^2,6*z2,2,0,0}};
coeffArr=LinearSolve[A,b];
gamma=x^Range[0,5].Reverse[coeffArr];
dgamma=D[gamma,x];
d2gamma=D[dgamma,x];
Print[format1@CForm[gamma]];
Print[format1@CForm[dgamma]];
Print[format1@CForm[d2gamma]];




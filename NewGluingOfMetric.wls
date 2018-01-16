#!/usr/bin/env wolframscript
(* ::Package:: *)

ClearAll[format1]
format1 = RawBoxes[Replace[ToBoxes@#, InterpretationBox[a_, b_, c___] :>
      With[{aa = StringReplace[a, 
            { "Sqrt" -> "sqrt", "Power(E," -> "exp(", "Power" -> "pow"}]}, aa], 
      {0, Infinity}]] &;
      (*Gluing of M(rho)*)
rho1=0.75;
tau1=1.0;
fp1=0.0;
fpp1=0.0;
rho2=1.25;
tau2=1.953125;
fp2=4.6875;
fpp2=7.5;
b={tau1,fp1,fpp1,tau2,fp2,fpp2};
A={{rho1^5,rho1^4,rho1^3,rho1^2,rho1,1},{5*rho1^4,4*rho1^3,3*rho1^2,2*rho1,1,0},{20*rho1^3,12*rho1^2,6*rho1,2,0,0},{rho2^5,rho2^4,rho2^3,rho2^2,rho2,1},{5*rho2^4,4*rho2^3,3*rho2^2,2*rho2,1,0},{20*rho2^3,12*rho2^2,6*rho2,2,0,0}};
coeffArr1=LinearSolve[A,b];
M=rho^Range[0,5].Reverse[coeffArr1];
Mp=D[M,rho];
Mpp=D[Mp,rho];
(*Gluing of t0(rho)*)
tau1=0.75;
fp1=1.0;
fpp1=0.0;
tau2=1.;
fp2=0.0;
fpp2=0.0;
b={tau1,fp1,fpp1,tau2,fp2,fpp2};
A={{rho1^5,rho1^4,rho1^3,rho1^2,rho1,1},{5*rho1^4,4*rho1^3,3*rho1^2,2*rho1,1,0},{20*rho1^3,12*rho1^2,6*rho1,2,0,0},{rho2^5,rho2^4,rho2^3,rho2^2,rho2,1},{5*rho2^4,4*rho2^3,3*rho2^2,2*rho2,1,0},{20*rho2^3,12*rho2^2,6*rho2,2,0,0}};
coeffArr2=LinearSolve[A,b];
t0=rho^Range[0,5].Reverse[coeffArr2];
t0p=D[t0,rho];
t0pp=D[t0p,rho];
X=(Mp*(t0-tau)+2*M*t0p)/((6*M^2*(t0-tau))^(1/3));
X=Simplify[X];
Y=(4.5*M*(t0-tau)^2)^(1/3);
Y=Simplify[Y];
dXdrho=D[X,rho];
dXdtau=D[X,tau];
dYdrho=D[Y,rho];
dYdtau=D[Y,tau];
dXdrho=Simplify[dXdrho];
dXdtau=Simplify[dXdtau];
dYdrho=Simplify[dYdrho];
dYdtau=Simplify[dYdtau];
Print["X:"];
Print[format1@CForm[X]];
Print["dX/drho"];
Print[format1@CForm[dXdrho]];
Print["dX/dtau"];
Print[format1@CForm[dXdtau]];
Print["Y"];
Print[format1@CForm[Y]];
Print["dY/drho"];
Print[format1@CForm[dYdrho]];
Print["dY/dtau"];
Print[format1@CForm[dYdtau]];











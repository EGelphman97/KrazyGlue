#!/usr/bin/env wolframscript
(* ::Package:: *)

z1=0.4;
w1=0.0;
fp1=0.0;
fpp1=0.0;
z2=1.7532659415312257;
w2=-0.8489326119968736;
fp2=-0.9304135317539163;
fpp2=-0.015209678948106337;
b={w1,fp1,fpp1,w2,fp2,fpp2};
A={{z1^5,z1^4,z1^3,z1^2,z1,1},{5*z1^4,4*z1^3,3*z1^2,2*z1,1,0},{20*z1^3,12*z1^2,6*z1,2,0,0},{z2^5,z2^4,z2^3,z2^2,z2,1},{5*z2^4,4*z2^3,3*z2^2,2*z2,1,0},{20*z2^3,12*z2^2,6*z2,2,0,0}};
coeffArr=LinearSolve[A,b];
gamma=x^Range[0,5].Reverse[coeffArr];
dgamma=D[gamma,x];
Print[gamma];
Print[dgamma];
Print[D[dgamma,x]];




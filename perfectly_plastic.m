% CODE TO FIND OPTIMUM AUTOFRETTAGE PRESSURE FOR PERFECTLY PLASTIC CONDITION
%MATERIAL PROPERTIES E=200GPA YIELD STRESS=250MPA
clear all
clc
close all
%Initiate dimensionless Elasto-Plastic boundary as a symbolic variable
syms rhoc
%Define Beta = Outer radius(=150mm)/Inner Radius(=60mm)
Beta=2.5;
%Define dimensionless Working pressure as = 250Mpa / Yield Stress(=250MPa)
Pw=1;
%Pressure applied as a function of rhoc
P=2*log(rhoc)+(2/3)*(1-(rhoc/Beta)^3);
%To Find Shakedown Pressure
%We need to evaluate parameters at inner radius rho=1
rho=1;
%Forward loading Hoop stress expression at inner radius
Stheta_forward = 2*log(rho)+1-P;
% Elastic unloading Hoop Stress expression at inner radius
Stheta_reverse=-P*(2*(rho)^3+(Beta)^3)/(2*rho^3*((Beta)^3-1));
%Residual Hoop Stress at inner radius
Stheta_residual = Stheta_forward+Stheta_reverse;
%Applying yield condition in reverse loading with residual radial stress
%zero at inner radius
yield_exp = Stheta_residual+1;
%Calculate maximum rhoc
rhocmax=vpasolve(yield_exp,1)
%Calculate Shakedown Pressure
P_max=(2*log(rhocmax)+(2/3)*(1-(rhocmax/Beta)^3))*250
%To Optimise Autofrettage pressure to be applied by method of total hoop
%stress optimisation
%Total hoop stress for a given working pressure at rhoc
Stheta_total=(2*log(rhoc)+1-P)+(2*rhoc^3+Beta^3)*(Pw-P)/(2*rhoc^3*(Beta^3-1));
%To optimise hoop stress to get optimum rhoc
exp=diff(Stheta_total);
%Calculate optimum rhoc
rhocopt=vpasolve(exp,1)
%Calculate optimum Autofrettage Pressure
Popt=(2*log(rhocopt)+(2/3)*(1-(rhocopt/Beta)^3))*250
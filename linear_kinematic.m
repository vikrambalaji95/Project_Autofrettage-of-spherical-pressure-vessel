clc
clear all
syms rhoc %Elasto-plastic boundary
m=0.02; %Ratio of Tangent Modulus and Elastic Modulus
mu = 0.3; %poisson's ratio
Beta = 2.5;%Ratio of outer radius to inner radius
Betac = Beta/rhoc;
Pw = 1; %Working pressure =250MPa
Denm = 2*m*(1-mu)+(1-m);
P1=(4/3)*(1-mu)*m*(1-(1/Beta^3))*rhoc^3;
P2=2*(1-m)*log(rhoc);
P3=(2/3)*(1-m)*(1-(1/Betac^3));
P=(P1+P2+P3)/Denm; %Pressure applied as a function of rhoc
rho=1; %at inner radius
%forward hoop stress at inner radius
Stheta_forward= -P+(2*(1-m)*log(rho)/Denm)+(((2/3)*(1-mu)*(rhoc^3/rho^3)*m*(2*rho^3+1))/Denm);
%forward radial stress at inner radius
Sradial_forward=-P+(2*(1-m)*log(rho)/Denm)+(((4/3)*(1-mu)*(rhoc^3/rho^3)*m*(rho^3-1))/Denm);
%forward equivalent stress at inner radius
S_forward=Stheta_forward-Sradial_forward;
%reverse hoop stress at inner radius
Stheta_revers=-1*P*(2*rho^3+Beta^3)/(2*rho^3*(Beta^3-1));
%residual hoop stress at inner radius
Stheta_residue=Stheta_forward+Stheta_revers;
expression=simplify(Stheta_residue);
%yield condition at onset of yielding at inner radius while unloading considering kinematic hardening
yield_exp=Stheta_residue+2-S_forward;
%function to solve rhoc at shakedown
rhocmax = vpasolve(yield_exp,1)
Betacmax=Beta/rhocmax;
P1_max=(4/3)*(1-mu)*m*(1-(1/Beta^3))*rhocmax^3;
P2_max=2*(1-m)*log(rhocmax);
P3_max=(2/3)*(1-m)*(1-(1/Betacmax^3));
P_max = 250*(P1_max+P2_max+P3_max)/Denm
%CODE FOR OBTAINING THE DIFFERENTIATED EQUATION AND TO DETERMINE OPTIMUM ROHC
% Forward(loading)hoop stress at rhoc
Stheta1= -P+(2*(1-m)*log(rhoc)/Denm)+(((2/3)*(1-mu)*m*(2*rhoc^3+1))/Denm);
% reverse(unloading)hoop stress at rhoc
Stheta_revers1=-P*(2*rhoc^3+Beta^3)/(2*rhoc^3*(Beta^3-1));
% reloading hoop stress at rhoc applying working pressure
Stheta_workload=Pw*(2*rhoc^3+Beta^3)/(2*rhoc^3*(Beta^3-1));
% total working hoop stress at rhoc
Stheta_total=Stheta1+Stheta_revers1+Stheta_workload;
expression1=simplify(Stheta_total);
%To optimise total working hoop stress to get optimum rhoc
expression2=diff(expression1);
%optimum rhoc
rhoc_opt = vpasolve(expression2,1)
Betac_opt=Beta/rhoc_opt;
P1_opt = (4/3)*(1-mu)*m*(1-(1/Beta^3))*rhoc_opt^3;
P2_opt = 2*(1-m)*log(rhoc_opt);
P3_opt=(2/3)*(1-m)*(1-(1/Betac_opt^3));
%optimum autofrettage pressure
P_opt = 250*(P1_opt+P2_opt+P3_opt)/Denm
clear

param.lambda = 0.87;
param.phi = 0.75;
param.w_m = 7;
param.w_f = param.phi*param.w_m;
param.gamma = 1.9;
param.cbar = 0.6;

A = [];
b = []; 
Aeq = [];
beq = [];


[l_star,v] = ...
    fmincon(@(l_f)-u_married_hh(param,l_f),0.5,A,b,Aeq,beq,0,1);

v = -v;
function [cost,l_star_init,l_star_final] = cost_function_test(data,param,x)

param.lambda = x(1);
param.cbar = x(2);

A = [];
b = []; 
Aeq = [];
beq = [];

[l_star_init,~] = ...
    fmincon(@(l_f)-u_married_hh(param,l_f),0.5,A,b,Aeq,beq,0,1);

param.w_m = data.growth;
param.w_f = param.phi*data.growth;

[l_star_final,~] = ...
    fmincon(@(l_f)-u_married_hh(param,l_f),0.5,A,b,Aeq,beq,0,1);

cost = 100*((l_star_init-data.l_init)^2+(l_star_final-data.l_final)^2);
%{
Computes the loss function to be minimized to select the weights in the
utility function for single households
%}

function L = timealloclossfn_s(param,targ,wage,x)

param.sigma_c = x(1);
param.sigma_l = x(2);
param.sigma_g = x(3);

[~,l,~,h,~,n] = per_sol_singles(param,wage);

L = (l-targ.l)^2+(h-targ.h)^2+(n-targ.n)^2;


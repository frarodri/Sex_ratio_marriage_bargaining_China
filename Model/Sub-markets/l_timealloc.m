%{
Computes the loss function that will be minimized to select the
parameters in the married household's utility function
%}

function L = l_timealloc(param,targ,wage_f,wage_m,pwf,x)

param.sigma_c= x(1);
param.sigma_l= x(2);
param.sigma_g = x(3);

[~,~,lf,lm,~,hf,hm,~,~,nf,nm] = ...
    per_sol_married(param,wage_f,wage_m,pwf,'quietly','true');
           
L = (lf-targ.lf)^2+(lm-targ.lm)^2+ (hf-targ.hf)^2+...
    (hm-targ.hm)^2+(nf-targ.nf)^2+(nm-targ.nm)^2;

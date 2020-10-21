%{
Computes the Pareto weight for the wife that implements the Egalitarian 
Bargaining solution within marriage, i.e. the Pareto weight that equalizes
the gains from marriage for both spouses, given a reservation quality
for marriage.
%}

function [Gf,Gm,pw] = egalbarg(param,wages,theta_S,qr,mu)

pw = ...
    fzero(@(pwf)mgains_diff(param,wages,theta_S,qr,pwf,mu),[10^-(10),1]);

[UMf,UMm] = uflow_married(param,wages.f,wages.m,pwf);

 
Gf = (UMf+cond_exp_q(qr,mu))/(1-param.beta*(1-param.delta));

Gm = (UMm+cond_exp_q(qr,mu))/(1-param.beta*(1-param.delta));


        
        
        
        
        
        





    
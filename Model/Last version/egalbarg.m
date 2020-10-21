%{
Computes the Pareto weight for the wife that implements the Egalitarian 
Bargaining solution within marriage, i.e. the Pareto weight that equalizes
the gains from marriage for both spouses, given a reservation quality
for marriage.
%}

function [Gf,Gm,pwf] = egalbarg(param,wage_f,wage_m,VSf,VSm,qr,mu)

pwf = fzero(@(pw)mgains_diff(param,wage_f,wage_m,VSf,VSm,qr,pw),...
    [10^-(10),1-10^-(10)]);

[UMf,UMm] = uflow_married(param,wage_f,wage_m,pwf,'quietly','true');

Gf = (1-normcdf(qr,mu))*(UMf+cond_exp_q(qr,mu))/...
    (1-param.beta*(1-param.delta));

Gm = (1-normcdf(qr,mu))*(UMm+cond_exp_q(qr,mu))/...
    (1-param.beta*(1-param.delta));


        
        
        
        
        
        





    
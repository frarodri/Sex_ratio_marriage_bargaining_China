%{
Given parameters, an expected sex ratio among singles, the flow utilities 
when married and single, a reservation match quality and a mean for the 
distribution of match quality draws, computes the value of being single for 
each sex
%}

function [VSf, VSm] = ...
    v_singles_submkts(param,theta_S,UMf,UMm,USf,USm,qr,mu)

[pi_f , pi_m] = match_prob(param,theta_S);

Nf = USf+param.beta*(1-param.delta)*pi_f*...
    ((1-normcdf(qr,mu))*(UMf+cond_exp_q(qr,mu))/...
    (1-param.beta*(1-param.delta)));

Df = 1-param.beta*(1-param.delta)*(1-pi_f*(1-normcdf(qr,mu)));

Nm = USm+param.beta*(1-param.delta)*pi_m*...
    ((1-normcdf(qr,mu))*(UMm+cond_exp_q(qr,mu))/...
    (1-param.beta*(1-param.delta)));

Dm = 1-param.beta*(1-param.delta)*(1-pi_m*(1-normcdf(qr,mu)));

VSf = Nf./Df;
VSm = Nm./Dm;
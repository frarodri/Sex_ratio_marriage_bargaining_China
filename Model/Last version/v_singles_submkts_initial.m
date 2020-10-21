%{
Given parameters, an expected sex ratio among singles, the flow utilities 
when married and single, a reservation match quality and a mean for the 
distribution of match quality draws, computes the value of being single for 
each sex
%}

function [VSf, VSm] = ...
    v_singles_submkts_initial(param,theta_S,USf,USm,UMf,UMm,qr,mu)

[pi_f , pi_m] = match_prob(param,theta_S);

PVSf = (USf+param.psi)/(1-param.beta*(1-param.delta));
PVSm = (USm)/(1-param.beta*(1-param.delta));

Nf = USf+param.psi+param.beta*(1-param.delta)*(param.rho*PVSf+...
    (1-param.rho)*pi_f*((1-normcdf(qr,mu))*(UMf+cond_exp_q(qr,mu))/...
    (1-param.beta*(1-param.delta))));

Df = 1-param.beta*(1-param.delta)*(1-param.rho)*...
    (1-pi_f*(1-normcdf(qr,mu)));

Nm = USm+param.beta*(1-param.delta)*(param.rho*PVSm+...
    (1-param.rho)*pi_m*((1-normcdf(qr,mu))*(UMm+cond_exp_q(qr,mu))/...
    (1-param.beta*(1-param.delta))));

Dm = 1-param.beta*(1-param.delta)*(1-param.rho)*...
    (1-pi_m*(1-normcdf(qr,mu)));

VSf = Nf./Df;
VSm = Nm./Dm;
%{
Given parameters, expectations, wages and reservation qualities for all
types of marriages (Q), computes the value of being single for each type
of agent and records it in two vectors, one for each sex
%}

function [VSf, VSm] = v_singles(param,expect,UMf,UMm,USf,USm,Q)

[pi_f , pi_m] = match_prob(param,expect.theta_S);

Nf = USf+param.beta*(1-param.delta)*pi_f*...
    ((1-normcdf(Q,param.MU)).*(UMf+cond_exp_q(param,Q))/...
    (1-param.beta*(1-param.delta)))*expect.Pm;

Df = 1-param.beta*(1-param.delta)*...
    (1-pi_f*(1-normcdf(Q,param.MU)*expect.Pm));

Nm = USm+param.beta*(1-param.delta)*pi_m*...
    ((1-normcdf(Q,param.MU)).*(UMm+cond_exp_q(param,Q))/...
    (1-param.beta*(1-param.delta)))'*expect.Pf;

Dm = 1-param.beta*(1-param.delta)*...
    (1-pi_m*(1-normcdf(Q,param.MU)'*expect.Pf));

VSf = Nf./Df;
VSm = Nm./Dm;
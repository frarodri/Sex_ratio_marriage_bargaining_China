%{
Computes the value of being single and record it in two vectors, one for 
each sex and each having as many entries as types there are for that sex
%}

function [VS_f, VS_m] = val_single(param,expect,wages,Q)

[pi_f , pi_m] = match_prob(param,expect.theta_S);

NUM_f = log(wages.OMEGA_f-param.cbar)+...
        param.beta*(1-param.delta)*pi_f*...
        ((1-normcdf(Q,param.MU)).*(flow_marr_chp(param,wages)+...
        cond_exp_q(param,Q))/(1-param.beta*(1-param.delta)))*expect.P_m;

NUM_m = log(wages.OMEGA_m-param.cbar)+...
        param.beta*(1-param.delta)*pi_m*...
        ((1-normcdf(Q,param.MU)).*(flow_marr_chp(param,wages)+...
        cond_exp_q(param,Q))/(1-param.beta*(1-param.delta)))'*expect.P_f;

DENOM_f = 1-param.beta*(1-param.delta)*...
    (1-pi_f*(1-normcdf(Q,param.MU)*expect.P_m));

DENOM_m = 1-param.beta*(1-param.delta)*...
    (1-pi_m*(1-normcdf(Q,param.MU)'*expect.P_f));

VS_f = NUM_f./DENOM_f;

VS_m = NUM_m./DENOM_m;
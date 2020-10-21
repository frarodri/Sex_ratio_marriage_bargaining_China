%{
Given parameters, expectations, wages and reservation qualities for all
types of marriages (Q), computes the value of being single for each type
of agent and records it in two vectors, one for each sex
%}

function [VSf, VSm] = v_singles(param,UMf,UMm,USf,USm,Q,theta_S)

types = size(param.Pf,1);
PIf = zeros(types,1);
PIm = zeros(types,1);

for z=1:types
   
    [pi_f , pi_m] = match_prob(param,theta_S(z));
    PIf(z) = pi_f;
    PIm(z) = pi_m;
    
end


Nf = USf+param.beta*(1-param.delta)*PIf.*...
    ((1-diag(normcdf(Q,param.MU))).*diag(UMf+cond_exp_q(Q,param.MU)))/...
    (1-param.beta*(1-param.delta));

Df = 1-param.beta*(1-param.delta)*...
    (1-PIf.*(1-diag(normcdf(Q,param.MU))));

Nm = USm+param.beta*(1-param.delta)*PIm.*...
    ((1-diag(normcdf(Q,param.MU))).*diag(UMm+cond_exp_q(Q,param.MU)))/...
    (1-param.beta*(1-param.delta));

Dm = 1-param.beta*(1-param.delta)*...
    (1-PIm.*(1-diag(normcdf(Q,param.MU))));

VSf = Nf./Df;
VSm = Nm./Dm;
function [theta_S_new,Q,PWf] = theta_update(param,wages,theta_S_prior)

%% Solve the marriage sub-markets
types = size(param.Pf,1);
Q_submkts = zeros(types,1);
PW_submkts = zeros(types,1);
VSf = zeros(types,1);
VSm = zeros(types,1);
PIf = zeros(types,1);
PIm = zeros(types,1);

for z=1:types
    
    mu = param.MU(z,z);
    theta_S = theta_S_prior(z);
    wage_f = wages.f(z);
    wage_m = wages.m(z);
    [q,pw,VSf_submkt,VSm_submkt] = ...
        submkt_pree(param,theta_S,wage_f,wage_m,mu,'quietly','true');
    [pi_f,pi_m] = match_prob(param,theta_S);
    
    Q_submkts(z) = q;
    PW_submkts(z) = pw;
    VSf(z) = VSf_submkt;
    VSm(z) = VSm_submkt;
    PIf(z) = pi_f;
    PIm(z) = pi_m;
    
end
%% Solve the entry markets

[Q,PWf] = entry_pree(param,wages,VSf,VSm,'quietly','true');

%% Compute the realized steady-state tightness for each sub-market

% Sf_0 = param.Pf.*(normcdf(Q,param.MU)*param.Pm);
% Sm_0 = param.theta_0*param.Pm.*...
%        ((param.theta_0-1)/param.theta_0+...
%        (1/param.theta_0)*(normcdf(Q,param.MU))'*param.Pf);
% theta_0 = Sm_0./Sf_0;
% 
% theta_S_new = theta_S_SS(param,theta_0,Q,theta_S_prior);

MRf = PIf.*(1-normcdf(Q_submkts,diag(param.MU)));
MRm = PIm.*(1-normcdf(Q_submkts,diag(param.MU)));

Sm = param.theta_0*param.Pm.*...
    ((param.theta_0-1)/param.theta_0+...
    (1/param.theta_0)*(normcdf(Q,param.MU))'*param.Pf)./...
    (1-(1-param.delta)*(1-MRm));
Sf = param.Pf.*(normcdf(Q,param.MU)*param.Pm)./...
    (1-(1-param.delta)*(1-MRf));

theta_S_new = Sm./Sf;

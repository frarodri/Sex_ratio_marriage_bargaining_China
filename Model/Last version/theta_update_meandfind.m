%{
Given a set of parameters, wages, targets and an initial vector of expected
sex ratios in each of the same type submarkets, computes the submarkets and
entry markets partial rational expectations equilibria. With these, it
calculates the means of the match quality draws among entrants to
replicate as closely as possible the marital sorting as measured by the
targeted contingency matrix. Finally, computes the steady state sex ratios
under the current expectations.
%}

function [theta_S_new,Q,PWf,MU] = ...
    theta_update_meandfind(param,wages,targ,theta_S_prior)

types = size(param.Pf,1);

% Solve the marriage sub-markets
Q_submkts = zeros(types,1);
PW_submkts = zeros(types,1);
VSf = zeros(types,1);
VSm = zeros(types,1);
PIf = zeros(types,1);
PIm = zeros(types,1);

MU_diag = zeros(5,5);
% MU_diag(1,1) = 1;
% MU_diag(4,4) = 0.5;
% MU_diag(5,5) = 1;

for z=1:types
    
    mu = MU_diag(z,z);
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

% Solve the entry markets
[Q,PWf] = entry_pree(param,wages,VSf,VSm,'quietly','true');

% Calculate the means for the draws that would most closely replicate the
% observed sorting
types = size(Q,1);
x0 = zeros(types*(types-1),1);

options = optimoptions('fminunc','Display','off');
mus = fminunc(@(MU_vec)l_maritalsort(param,targ,theta_S_prior,Q,MU_vec),...
    x0,options);

MU_aux = vec2mat(mus,types);
MU_l = [zeros(1,types);tril(MU_aux)];
MU_u = [triu(MU_aux,1);zeros(1,types)];
MU = MU_l+MU_u+MU_diag;
param.MU = MU;

% Compute steady-state tightness for each sub-market
MRf = PIf.*(1-normcdf(Q_submkts,diag(param.MU)));
MRm = PIm.*(1-normcdf(Q_submkts,diag(param.MU)));

Sm = param.theta_0*param.Pm.*...
    ((param.theta_0-1)/param.theta_0+(normcdf(Q,param.MU))'*param.Pf)./...
    (1-(1-param.delta)*(1-MRm));
Sf = param.Pf.*(normcdf(Q,param.MU)*param.Pm)./...
    (1-(1-param.delta)*(1-MRf));

theta_S_new = Sm./Sf; 

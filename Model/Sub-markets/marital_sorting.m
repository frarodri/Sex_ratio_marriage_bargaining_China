%{
Given a set of parameters (including notably the means of the distribution
for the match quality draws, a vector of sex ratios for each of the 
same-type marriage markets, and a matrix of reservation match qualities for 
marriage, computes the marital sorting matrix
%}

function MS = marital_sorting(param,theta_S,Q)

types=size(param.Pf,1);

% Calculate the total measure of each combination-type of marriages that
% takes place between entrants
M_x = param.Pf*param.Pm'.*(1-normcdf(Q,param.MU))/param.delta;

% Compute the marriage rate in the same-type markets 
PIf = zeros(types,1);
PIm = zeros(types,1);

for z=1:types
    
    theta_Sz = theta_S(z);
    [pi_f , pi_m] = match_prob(param,theta_Sz);
    PIf(z) = pi_f;
    PIm(z) = pi_m;
    
end

MRf = PIf.*(1-normcdf(diag(Q),diag(param.MU)));
MRm = PIm.*(1-normcdf(diag(Q),diag(param.MU)));

% Calculate the measure of singles in the same-type marriage markets
Sm = param.theta_0*param.Pm.*...
    ((param.theta_0-1)/param.theta_0+(normcdf(Q,param.MU))'*param.Pf)./...
    (1-(1-param.delta)*(1-MRm));
Sf = param.Pf.*(normcdf(Q,param.MU)*param.Pm)./...
    (1-(1-param.delta)*(1-MRf));

% Compute the measure of same-type marriages that happen after entry
% (within same-type markets)
M_zz_f = Sf.*MRf/param.delta;
M_zz_m = Sm.*MRm/param.delta;
M_zz = (M_zz_m+M_zz_f)/2;

% Compute matrix with total measure of each type of marriage
MS = M_x+diag(M_zz);

% Compute marital sorting matrix
MS = MS/sum(sum(MS));











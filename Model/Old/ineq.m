function [ineq_inc, ineq_cons] = ineq(param,wages,Q,MS,theta_S,P_f,P_m)

%{
Compute the Gini index for household income and consumption
%}

% Measure of married couples
[pi_f, ~] = match_prob(param,theta_S);

mr_f = pi_f*P_f'*(1-normcdf(Q,param.MU))*P_m;
%mr_m = pi_m*Z_f'*(1-normcdf(Q,param.MU))*Z_m;

married_measure = mr_f/param.delta;

married_measure = reshape(married_measure*MS,[],1);

% Income of married couples
szOMEGA_f = size(wages.OMEGA_f,1);
szOMEGA_m = size(wages.OMEGA_m,1);

Lf = married_fls(param,wages);

married_inc = (ones(szOMEGA_f,szOMEGA_m).*...
    wages.OMEGA_m'+Lf.*wages.OMEGA_f)/2;

married_inc = reshape(married_inc,[],1);

% Consumption of married couples
married_cons = (ones(szOMEGA_f,szOMEGA_m).*...
    wages.OMEGA_m'+Lf.*wages.OMEGA_f)/1.7;

married_cons = reshape(married_cons,[],1);

% Measure of singles
singles_measure = [P_f ; theta_S*P_m];

% Income of singles
singles_inc = [wages.OMEGA_f ; wages.OMEGA_m];

% Concatenate measure vectors and income vectors
measures = [married_measure ; singles_measure];

% Concatenate income vectors
inc = [married_inc ; singles_inc];

% Concatenate consumption vectors
cons = [married_cons ; singles_inc];

% Compute the gini index
ineq_inc = gini(measures,inc);
ineq_cons = gini(measures,cons);





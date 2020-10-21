function MG = mgains_initial(param,wage_f,wage_m,USf,USm,theta_S,pwf,mu,qr)

% Compute the flow value of marriage with current Pareto weight for wife
[UMf,UMm] = uflow_married(param,wage_f,wage_m,pwf,'quietly','true');

% Compute the value of being single with current Pareto weight for wife
[VSf, VSm] = ...
    v_singles_submkts_initial(param,theta_S,UMf,UMm,USf,USm,qr,mu);

% Compute the gains from marriage for wife and husband and the difference 
% with current Pareto weight for wife
Wf = (UMf+qr)/(1-param.beta*(1-param.delta))-VSf;
Wm = (UMm+qr)/(1-param.beta*(1-param.delta))-VSm;
MG = Wf+Wm;
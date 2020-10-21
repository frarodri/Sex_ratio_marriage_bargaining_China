function NP = nash_product(param,wages,theta_S,qr,pwf,mu)

% Compute the flow value of being single
[USf,USm] = uflow_singles(param,wages);

% Compute the flow value of marriage with current Pareto weight for wife
[UMf,UMm] = uflow_married(param,wages.f,wages.m,pwf,'quietly','true');

% Compute the value of being single with current Pareto weight for wife
[VSf, VSm] = v_singles_submkts(param,theta_S,UMf,UMm,USf,USm,qr,mu);

% Compute the gains from marriage for wife and husband 
Wf = (UMf+qr)/(1-param.beta*(1-param.delta))-VSf;
Wm = (UMm+qr)/(1-param.beta*(1-param.delta))-VSm;

% Compute the Nash product
NP = Wf*Wm;


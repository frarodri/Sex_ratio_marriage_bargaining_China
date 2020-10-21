function NP = nash_product_entrants(param,wages,VSf,VSm,qr,pwf)

% Compute the flow value of marriage with current Pareto weight for wife
[UMf,UMm] = uflow_married(param,wages.f,wages.m,pwf,'quietly','true');

% Compute the gains from marriage for wife and husband 
Wf = (UMf+qr)/(1-param.beta*(1-param.delta))-VSf;
Wm = (UMm+qr)/(1-param.beta*(1-param.delta))-VSm;

% Compute the Nash product
NP = Wf*Wm;


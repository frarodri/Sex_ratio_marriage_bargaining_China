function DW = mgains_diff(param,wage_f,wage_m,VSf,VSm,qr,pwf)

% Compute the flow value of marriage with current Pareto weight for wife
[UMf,UMm] = uflow_married(param,wage_f,wage_m,pwf,'quietly','true');

% Compute the gains from marriage for wife and husband and the difference 
% with current Pareto weight for wife
Wf = (UMf+qr)/(1-param.beta*(1-param.delta))-VSf;
Wm = (UMm+qr)/(1-param.beta*(1-param.delta))-VSm;
DW = Wf-Wm;

function DW = mgains_diff(param,expect,wages,q,pw)

% Compute the flow value of being single
[USf,USm] = uflow_singles(param,wages);

% Compute the flow value of marriage with current Pareto weight for wife
[UMf,UMm] = uflow_married(param,wages,pw,'quietly','true');

% Compute the value of being single with current Pareto weight for wife
[VSf, VSm] = v_singles(param,expect,UMf,UMm,USf,USm,q);

% Compute the gains from marriage and the difference with current Pareto 
% weight for wife
Wf = (UMf+q)/(1-param.beta*(1-param.delta))-VSf;
Wm = (UMm+q)/(1-param.beta*(1-param.delta))-VSm;
DW = Wf-Wm;
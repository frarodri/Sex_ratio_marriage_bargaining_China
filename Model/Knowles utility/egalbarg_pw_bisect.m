%{
Computes the Pareto weight for the wife that implements the Egalitarian 
Bargaining solution within marriage, i.e. the Pareto weight that equalizes
the gains from marriage for both spouses, given a reservation quality
for marriage. The algortihm used is the bisection method.
%}

function mu_f = egalbarg_pw_bisect(param,theta,wage_f,wage_m,q)

% Set wages and expectations
wages.f = wage_f;
wages.m = wage_m;

expect.Pf = 1;
expect.Pm = 1;
expect.theta_S = theta;

% Initialize wife's Pareto weight to 0.5 and upper and lower bounds to 
% 1 and 0, respectively 
mu_f = 0.5;
ub = 1;
lb = 0;

% Compute the flow value of being single
[USf,USm] = uflow_singles(param,wages);

% Compute the flow value of marriage with current Pareto weight for wife
[UMf,UMm] = uflow_married(param,wages,mu_f);

% Compute the value of being single with current Pareto weight for wife
[VSf, VSm] = v_singles(param,expect,UMf,UMm,USf,USm,q);

% Compute the gains from marriage and the difference with current Pareto 
% weight for wife
Wf = (UMf+q)/(1-param.beta*(1-param.delta))-VSf;
Wm = (UMm+q)/(1-param.beta*(1-param.delta))-VSm;
DW = Wf-Wm;

% Set tolerance level for convergence
tol = 10^(-10);

% Start bisection algortihm

while abs(DW)>tol
    
    if DW>0
        
        ub = mu_f;
        mu_f = (0.5+lb)/2;
                
    else
        
        lb = mu_f;
        mu_f = (0.5+ub)/2;
                
    end
               
    % Re-compute the flow value of marriage with current Pareto weight
    % for wife
    [UMf,UMm] = uflow_married(param,wages,mu_f);


    % Re-compute the value of being single with current Pareto weight 
    % for wife
    [VSf, VSm] = v_singles(param,expect,UMf,UMm,USf,USm,q);

    % Re-compute the gains from marriage and the difference with 
    % current Pareto weight for wife
    Wf = (UMf+q)/(1-param.beta*(1-param.delta))-VSf;
    Wm = (UMm+q)/(1-param.beta*(1-param.delta))-VSm;
    DW = Wf-Wm;
        
end
        
        
        
        
        
        





    
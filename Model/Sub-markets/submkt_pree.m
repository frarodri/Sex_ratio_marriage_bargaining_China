%{
Computes the Partial Rational Expectations Equilibrium (PREE) with
Egalitarian Bargaining for a marriage sub-market, i.e. a market in which 
all agents have the same type. That is, given a set of parameters, an 
expected sex ratio among singles, and wages for females and males, computes 
the reservation match quality and the Pareto weight of the wife under 
Egalitarian Bargaining for marriage between agents with the same type. 
Computes also the value of being single for each sex, which serves as a 
continuation value for the problem of entrants that are presented with a 
marriage opportunity with an agent that potentially has a different type.
%}

function [qr,pw,VSf,VSm] = ...
    submkt_pree(param,theta_S,wage_f,wage_m,mu,varargin)

p = inputParser;
paramName = 'quietly';
defaultVal = "false";
errorMsg = 'Quietly must be either true or false';
validationFcn = @(x) assert(strcmp(x,"true") | strcmp(x,"false"),errorMsg);
addParameter(p,paramName,defaultVal,validationFcn);
parse(p,varargin{:});

% Guess an initial reservation match qualitiy q
qr = 0;

% Set wages 
wages.f = wage_f;
wages.m = wage_m;

% Compute the flow value of being single
[USf,USm] = uflow_singles(param,wages);
 
% Set initial value for error and tolerance for loop
current_error = 1;
tol = 10^(-5);
it = 1;

while current_error>tol
   
    qr_prior = qr;
    
    % 1. Compute the Pareto weights that implement the Egalitarian 
    % Bargaining solution for the current reservation match quality
%     pw = egalbarg_pw(param,theta_S,wage_f,wage_m,qr,mu);
    pw = nashbarg_pw(param,theta_S,wage_f,wage_m,qr,mu);
    
    % 2. Compute the flow value of marriage for the agents with the current
    % Pareto weights
    [UMf,UMm] = uflow_married(param,wage_f,wage_m,pw,'quietly','true');
    
    % 3. Compute the value of being single for the agents with the current
    % Pareto weights and reservation match quality
    [VSf,VSm] = v_singles_submkts(param,theta_S,UMf,UMm,USf,USm,qr,mu);
    
    % 4. Update the reservation match quality to set the gains from
    % marriage to zero
%     qf = (1-param.beta*(1-param.delta))*VSf-UMf;
%     qm = (1-param.beta*(1-param.delta))*VSm-UMm; % Old method
%     qr = (qf+qm)/2;
%     error_pw = abs(qf-qm);
    qr = fzero(@(qr)mgains(param,wages,theta_S,qr,pw,mu),0);
    error_report = abs(qr-qr_prior);
    
    if strcmp(p.Results.quietly,"false")
               
%         fprintf(['Iteration %d, Pareto weight error is %g\n, current '...
%             'error is %g\n'],it,error_pw,error_report);

        fprintf(['Iteration %d, current error is %g\n'],it,error_report);

    end
    
    current_error = error_report;
    it = it+1;
    
end
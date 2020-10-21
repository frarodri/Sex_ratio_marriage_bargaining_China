%{
Computes the Partial Rational Expectations Equilibrium (PREE) with
bargaining solution for a marriage market in which all agents have the same
type. That is, given a set of parameters, an expected sex ratio among
singles, wages, and the mean of the distribution of the match quality
draws, computes the reservation match quality for marriage and the Pareto
weight for the wife that implements the selected bargaining solution. Also
computes the value of being single, which serves as a continuation value in
the problem faced by new entrants.
%}
  
function [qr,pwf,VSf,VSm] = ...
    submkt_pree(param,theta_S,wage_f,wage_m,mu,USf,USm,start_pw,varargin)

p = inputParser;
paramName1 = 'Display';
defaultVal1 = 'all';
errorMsg1 = 'Display must be iter, time, all or quiet';
validationFcn1 = ...
    @(x) assert(strcmp(x,'iter') | strcmp(x,'time') | strcmp(x,'all') ...
    | strcmp(x,'quiet'),errorMsg1);
addParameter(p,paramName1,defaultVal1,validationFcn1);

paramName2 = 'bargaining';
defaultVal2 = "Nash";
errorMsg2 = 'The bargaining solution must be Egalitarian or Nash';
validationFcn2 = ...
    @(x) assert(strcmp(x,"Egalitarian") | strcmp(x,"Nash"),errorMsg2);
addParameter(p,paramName2,defaultVal2,validationFcn2);

parse(p,varargin{:});

%% Preliminaries
% Set initial Pareto weight of wife and compute the initial indirect
% utility flow under this Pareto weight for married people
pwf = start_pw;
[UMf,UMm] = uflow_married(param,wage_f,wage_m,pwf,'quietly','true');

% Compute the initial reservation match quality and value of being single
qr = fzero(...
    @(q)mgains_initial(param,wage_f,wage_m,USf,USm,theta_S,pwf,mu,q),mu);

[VSf, VSm] = ...
    v_singles_submkts_initial(param,theta_S,USf,USm,UMf,UMm,qr,mu);

% Set initial value for error and tolerance for loop
current_error = 1;
tol = 10^(-5);
it = 1;

%% Loop that computes the PREE
t = cputime;
while current_error>tol
 
    qr_prior = qr;
    pwf_prior = pwf;
    
    % 1. Update the Pareto weights and the value of being single according
    % to the bargaining rule
    
    if strcmp(p.Results.bargaining,"Nash")
        
        [Gf,Gm,pwf] = ...
            nashbarg(param,wage_f,wage_m,VSf,VSm,qr,mu,pwf_prior);
        [VSf, VSm] = ...
            v_singles_submkts(param,theta_S,USf,USm,Gf,Gm,qr,mu);
        
    else
       
        [Gf,Gm,pwf] = egalbarg(param,wage_f,wage_m,VSf,VSm,qr,mu);
        [VSf, VSm] = ...
            v_singles_submkts(param,theta_S,USf,USm,Gf,Gm,qr,mu);
        
    end
    
    % 2. Update the reservation match quality by setting the gains from
    % marriage to zero
    
    qr = fzero(@(q)mgains(param,wage_f,wage_m,VSf,VSm,pwf,q),qr_prior);
    
    % 3. Update the error
    
    error_report = max(abs(qr-qr_prior),abs(pwf-pwf_prior));
    
    if strcmp(p.Results.Display,'iter') || strcmp(p.Results.Display,'all')
        
        message = ['Iteration %d, Pareto weight of wife was %g,'...
            'reservation match quality was %g, current error is %g\n'];
        fprintf(message,it,pwf,qr,error_report);
         
    end
    
    current_error = error_report;
    it = it+1;
    
end
e = cputime-t;

if strcmp(p.Results.Display,'time') || strcmp(p.Results.Display,'all')

     fprintf('Submarket loop elapsed time: %g seconds\n',e);

end


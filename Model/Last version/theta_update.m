function [theta_S_new,Q,PWf] = ...
    theta_update(param,wages,theta_S_prior,USf,USm,start_PWf,varargin)

p = inputParser;
paramName1 = 'quietly';
defaultVal1 = "false";
errorMsg1 = 'Quietly must be either true or false';
validationFcn1 = ...
    @(x) assert(strcmp(x,"true") | strcmp(x,"false"),errorMsg1);
addParameter(p,paramName1,defaultVal1,validationFcn1);

paramName2 = 'bargaining';
defaultVal2 = "Nash";
errorMsg2 = 'The bargaining solution must be Egalitarian or Nash';
validationFcn2 = ...
    @(x) assert(strcmp(x,"Egalitarian") | strcmp(x,"Nash"),errorMsg2);
addParameter(p,paramName2,defaultVal2,validationFcn2);

paramName3 = 'Display';
defaultVal3 = 'all';
errorMsg3 = 'Display must be iter, time, all or quiet';
validationFcn3 = ...
    @(x) assert(strcmp(x,'iter') | strcmp(x,'time') | strcmp(x,'all') ...
    | strcmp(x,'quiet'),errorMsg3);
addParameter(p,paramName3,defaultVal3,validationFcn3);

parse(p,varargin{:});

t = cputime;

%% Solve the marriage sub-markets
types = size(param.Pf,1);
Q_submkts = zeros(types,1);
PW_submkts = zeros(types,1);
VSf = zeros(types,1);
VSm = zeros(types,1);
PIf = zeros(types,1);
PIm = zeros(types,1);

for z=1:types
    
    mu = param.MU(z,z);
    theta_S = theta_S_prior(z);
    wage_f = wages.f(z);
    wage_m = wages.m(z);
    us_f = USf(z);
    us_m = USm(z);
    start_pw = start_PWf(z);
    
    [q,pw,VSf_submkt,VSm_submkt] = ...
        submkt_pree(param,theta_S,wage_f,wage_m,mu,us_f,us_m,start_pw,...
        'Display',p.Results.Display,'bargaining',p.Results.bargaining);
    
    [pi_f,pi_m] = match_prob(param,theta_S);
    
    Q_submkts(z) = q;
    PW_submkts(z) = pw;
    VSf(z) = VSf_submkt;
    VSm(z) = VSm_submkt;
    PIf(z) = pi_f;
    PIm(z) = pi_m;
    
end
%% Solve the entry markets

[Q,PWf] = entry_pree(param,wages,VSf,VSm,...
    'quietly',p.Results.quietly,'bargaining',p.Results.bargaining);

%% Compute the realized steady-state tightness for each sub-market

% Sf_0 = param.Pf.*(normcdf(Q,param.MU)*param.Pm);
% Sm_0 = param.theta_0*param.Pm.*...
%        ((param.theta_0-1)/param.theta_0+...
%        (1/param.theta_0)*(normcdf(Q,param.MU))'*param.Pf);
% theta_0 = Sm_0./Sf_0;
% 
% theta_S_new = theta_S_SS(param,theta_0,Q,theta_S_prior);

MRf = PIf.*(1-normcdf(Q_submkts,diag(param.MU)));
MRm = PIm.*(1-normcdf(Q_submkts,diag(param.MU)));

Sm = param.theta_0*param.Pm.*...
    ((param.theta_0-1)/param.theta_0+...
    (1/param.theta_0)*(normcdf(Q,param.MU))'*param.Pf)./...
    (1-(1-param.delta)*(1-param.rho)*(1-MRm));
Sf = param.Pf.*(normcdf(Q,param.MU)*param.Pm)./...
    (1-(1-param.delta)*(1-param.rho)*(1-MRf));

theta_S_new = Sm./Sf;

e = cputime-t;

if strcmp(p.Results.quietly,"false")

    message_elapsed_time = ['Elapsed time for the update of the sex '...
        'ratio among singles: %g seconds \n'];
    fprintf(message_elapsed_time,e);

end

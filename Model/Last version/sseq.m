function [theta_S,Q,PWf] = sseq(param,wages,USf,USm,varargin)

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

current_err = 1;
tol = 10^(-3);
it = 1;

theta_S_prior=param.Pm*param.theta_0./param.Pf;

types = size(param.Pf,1);
start_PWf = (1/2)*ones(types,1);

t = cputime;
while current_err>tol && it<=50
    
    [theta_S,Q,PWf] = theta_update(param,wages,theta_S_prior,USf,USm,...
        start_PWf,...
        'display',p.Results.Display,...
        'quietly',p.Results.quietly,...
        'bargaining',p.Results.bargaining);
    error_report = sum(abs(theta_S-theta_S_prior));
    theta_S_prior = theta_S;
    start_PWf = diag(PWf);
        
    if strcmp(p.Results.quietly,"false")
        
         fprintf('Iteration %d, current error is %g\n',it,error_report);
         
    end
    
    current_err = error_report;
    it = it+1;
    
end

e = cputime-t;

if strcmp(p.Results.quietly,"false")
        
    message_elapsed_time = ['Elapsed time until SS equilibrium was '...
        'found: %g seconds \n'];
    fprintf(message_elapsed_time,e);
    
end
    
    

    
    
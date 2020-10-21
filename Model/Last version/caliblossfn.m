%{
Computes the loss function to be minimized for the calibration
%}

function L = caliblossfn(param,wages,targ,USf,USm,calib_params,varargin)

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

parse(p,varargin{:});

%% Define the parameters to calibrate

types = size(param.Pf,1);

param.sigma_c = 1/exp(calib_params(1));
param.sigma_l = 1/exp(calib_params(2));
param.sigma_g = 1-param.sigma_c-param.sigma_l;

param.psi = calib_params(3);

param.MU = vec2mat(calib_params(4:end),types);

%% Compute the steady-state equilibrium

[theta_S,Q,PWf] = sseq(param,wages,USf,USm,...
    'quietly',p.Results.quietly,'bargaining',p.Results.bargaining);

%% Compute model moments

MS = marital_sorting(param,theta_S,Q);

[~,~,Lf,Lm,~,Hf,Hm,~,~,Nf,Nm] = ...
                      per_sol_married(param,wages.f,wages.m,PWf,'quietly','true');
                 
lf= sum(sum(MS.*Lf));
lm = sum(sum(MS.*Lm));
hf = sum(sum(MS.*Hf));
hm = sum(sum(MS.*Hm));
nf = sum(sum(MS.*Nf));
nm = sum(sum(MS.*Nm));

%% Compute the loss

L = (lf-targ.lf)^2+(lm-targ.lm)^2+(hf-targ.hf)^2+(hm-targ.hm)^2+...
    (nf-targ.nf)^2+(nm-targ.nm)^2+sum(sum((MS-targ.MS).^2));


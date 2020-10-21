function [UMf,UMm] = uflow_married(param,wages,PWf,varargin)

p = inputParser;
paramName = 'quietly';
defaultVal = "false";
errorMsg = 'Quietly must be either true or false';
validationFcn = @(x) assert(strcmp(x,"true") | strcmp(x,"false"),errorMsg);
addParameter(p,paramName,defaultVal,validationFcn);
parse(p,varargin{:});

if strcmp(p.Results.quietly,"true")
    
    [Cf,Cm,Lf,Lm,G] = per_sol_married(param,wages,PWf,'quietly','true');
    UMf = ind_u(param,Cf,Lf,G);
    UMm = ind_u(param,Cm,Lm,G);
    
else
    
    [Cf,Cm,Lf,Lm,G] = per_sol_married(param,wages,PWf);
    UMf = ind_u(param,Cf,Lf,G);
    UMm = ind_u(param,Cm,Lm,G);

end
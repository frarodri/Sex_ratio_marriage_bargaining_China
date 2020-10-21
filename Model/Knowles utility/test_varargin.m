function u = test_varargin(x1,x2,param,varargin)

p = inputParser;
paramName = 'quietly';
defaultVal = "false";
errorMsg = 'Quietly must be either true or false';
validationFcn = @(x) assert(strcmp(x,"true") | strcmp(x,"false"),errorMsg);
addParameter(p,paramName,defaultVal,validationFcn);
parse(p,varargin{:});

u = (x1.^(1-param.sigma))/(1-param.sigma)+...
    (x2.^(1-param.sigma))/(1-param.sigma);

if strcmp(p.Results.quietly,'false')
    
    disp('Done!');
    
end


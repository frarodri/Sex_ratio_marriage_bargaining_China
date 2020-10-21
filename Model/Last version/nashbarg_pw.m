%{
Computes the Pareto weight for the wife that implements the Nash 
Bargaining solution within marriage, i.e. the Pareto weight that maximizes 
the Nash product , given the values of being single. 
%}

function pw = nashbarg_pw(param,wage_f,wage_m,VSf,VSm,q,varargin)

p = inputParser;
paramName1 = 'start_point';
defaultVal1 = 0.5;
errorMsg1 = 'The start value must lie between 0 and 1';
validationFcn1 = ...
    @(x) assert(x>0 && x<1,errorMsg1);
addParameter(p,paramName1,defaultVal1,validationFcn1);

parse(p,varargin{:});

A = [];
b = []; 
Aeq = [];
beq = [];
x0 = rand();

if ~isempty(varargin)
    
    x0 = p.Results.start_point;
    
end

options = optimoptions('fmincon','Display','off');
pw = fmincon(@(pwf)-nash_product(param,wage_f,wage_m,VSf,VSm,q,pwf),...
                               x0,A,b,Aeq,beq,0,1,[],options);
                           
%{
Computes the partial rational expectations equilibrium (PREE) for the 
two-sided heterogeneous agents marriage seach model, i.e., given a set 
of parameters (param), expectations (expect) and wages (wages), computes 
a function q(zf,zm) of reservation match quality values for marriage 
between female agents of type zf and male agents of type zm. The main 
result is stored in matrix Q, in which rows represent female types and 
columns male types.
%}

function [Q,VSf,VSm,ID] = marriage_pree(param,expect,wages,PWf,varargin)

p = inputParser;
paramName = 'quietly';
defaultVal = "false";
errorMsg = 'Quietly must be either true or false';
validationFcn = @(x) assert(strcmp(x,"true") | strcmp(x,"false"),errorMsg);
addParameter(p,paramName,defaultVal,validationFcn);
parse(p,varargin{:});

szPf = size(expect.Pf,1);
szPm = size(expect.Pm,1);

szwagesf = size(wages.f,1);
szwagesm = size(wages.m,1);

dist_tol = 10^(-10);

if iscolumn(expect.Pf) && iscolumn(expect.Pm) && ...
   iscolumn(wages.f) && iscolumn(wages.m) && ...
   (abs(sum(expect.Pf)-1)<dist_tol) && (abs(sum(expect.Pm)-1)<dist_tol) && ...
   all(expect.Pf>=0) && all(expect.Pm>=0) && ...
   all(wages.f>0) && all(wages.m>0)
 
    if szPf == szwagesf && szPm == szwagesm
        
        % 1. Compute the utility flow of married people ex-companionship
        [UMf,UMm] = uflow_married(param,wages,PWf);
        
        % 2. Compute the utility flow of singles
        [USf,USm] = uflow_singles(param,wages);
                
        % 2. Guess an initial reservation match qualitiy matrix Q
        Q = zeros(szPf,szPm);
        
        % 3. Set initial value for error and tolerance for loop
        current_error = 1;
        tol = 10^(-10);
        it = 1;
        
        while current_error>tol
            
            Qpr = Q;
            
            % 4. Compute the value of being single for each sex and type
            % with the current reservation match qualities
            [VSf,VSm] = v_singles(param,expect,UMf,UMm,USf,USm,Q);
            
            % 5. Compute the value of marriage with the current reservation
            % qualitities
            VMf = (UMf+Q)/(1-param.beta*(1-param.delta));
            VMm = (UMm+Q)/(1-param.beta*(1-param.delta));
            
            % 6. Compute the surplus from marriage with the current 
            % reservation qualities
            Wf = VMf-VSf.*ones(szPf,szPm);
            Wm = VMm-VSm'.*ones(szPf,szPm);
            
            % 7. For each pair of female-male types, compute who has the 
            % smaller surplus from marriage under the current reservation
            % qualities
            ID = (Wf<=Wm);
            
            % 8. Compute a new Q by equalizing the value of marriage to the
            % value of staying single for each pair of female-male types.
            % Repeat 4-6 until convergence
            Q = ID.*((1-param.beta*(1-param.delta))*VSf.*ones(szPf,szPm)...
                -UMf)+(1-ID).*((1-param.beta*(1-param.delta))*...
                VSm'.*ones(szPf,szPm)-UMm);
            err_report = sum(sum(abs(Q-Qpr)));
            
            if strcmp(p.Results.quietly,"false")
                
                fprintf('Iteration %d, current error is %g\n',it,...
                    err_report);
                
            end
                
            current_error = err_report;
            it = it+1;
            
        end
        
    else
        
        message1 = ['Expected distribution of types and wage vectors '...
            'must have the same dimensions for each sex'];
        error(message1);
        
    end
    
else
    
    message2 = ['Expected distribution of types must be a column '...
        'vector with non-negative elements that sum to 1, and the '...
        'wages must be column vectors with positive entries'];
    error(message2);
    
end

end
%{
Computes (a) steady-state equilibrium under full rational expectations for 
the marriage market given a set of parameters (including the sex ratio and 
distributions over types among entrants), wages and a set of initial 
expectations  
%}

function [Q,theta_S,SPf,SPm] = ...
    marriage_sseq(param,wages,init_expect,PWf,varargin)

p = inputParser;
paramName = 'quietly';
defaultVal = "false";
errorMsg = 'Quietly must be either true or false';
validationFcn = @(x) assert(strcmp(x,"true") | strcmp(x,"false"),errorMsg);
addParameter(p,paramName,defaultVal,validationFcn);
parse(p,varargin{:});

szPf = size(param.Pf,1);
szPm = size(param.Pm,1);

szwagesf = size(wages.f,1);
szwagesm = size(wages.m,1);

if iscolumn(param.Pf) && iscolumn(param.Pm) && ...
   iscolumn(wages.f) && iscolumn(wages.m) && ...
   sum(param.Pf)==1 && sum(param.Pm)==1 && ...
   all(param.Pf>=0) && all(param.Pm>=0) && ...
   all(wages.f>0) && all(wages.m>0)
 
    if szPf == szwagesf && szPm == szwagesm
        
        expect.theta_S = init_expect.theta_S;
        expect.Pf = init_expect.Pf;
        expect.Pm = init_expect.Pm;
        
        current_error = 1; % Initialize error
        tol = 1*exp(-10); % Tolerance
        
        it = 1;    
        
        while current_error>tol
            
            % Compute marriage reservation quality matrix from PREE under
            % the current expectations
            [Q,~,~,~] = ...
                marriage_pree(param,expect,wages,PWf,'quietly','true');
            
            % Calculate the steady-state distributions resulting from the
            % current PREE  
            [theta_S, SPf, SPm] = dist_singles(param,Q);
            
            % Compute the new current error, i.e. the distance in absolute
            % value between the expected and steady-state sex ratio and 
            % distribution over types
            error_report = abs(theta_S-expect.theta_S)+...
                           sum(abs(SPf-expect.Pf))+...
                           sum(abs(SPm-expect.Pm));
                
            % Update the expectations to the steady-state equivalents
            expect.theta_S = theta_S;
            expect.Pf = SPf;
            expect.Pm = SPm;  
            
            if strcmp(p.Results.quietly,"false")
                
                fprintf('Iteration %d, current error is %g\n',it,...
                    error_report);
                
            end
            
            current_error = error_report;
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


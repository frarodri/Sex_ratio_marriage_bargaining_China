%{
Computes the partial rational expectations equilibrium (PREE) for the 
two-sided heterogeneous agents marriage seach model, i.e., given a set 
of parameters (param), expectations (expect) and wages (wages), computes 
a function q(z_f,z_m) of cut-off qualities for marriage between a female
agent of type z_f and a male agent of type z_m. The main result is stored
in matrix Q, in which rows represent female types and columns male types.
%}

function [Q, VS_f, VS_m, ID, MS] = marriage_pree(param,expect,wages)

szP_f = size(expect.P_f,1);
szP_m = size(expect.P_m,1);

szOMEGA_f = size(wages.OMEGA_f,1);
szOMEGA_m = size(wages.OMEGA_m,1);

if iscolumn(expect.P_f) && iscolumn(expect.P_m) && ...
        iscolumn(wages.OMEGA_m) && iscolumn(wages.OMEGA_m)
    
    if szP_f==szOMEGA_f && szP_m==szOMEGA_m
         
        % 1. Guess a matrix of cutoff match qualities, Q
        Q = zeros(szP_f,szP_m);

        err = 1; % Initialize error
        tol = 1*exp(-10); % Tolerance

        % 2. Compute the value of marriage ex companionship
        VM = flow_marr_chp(param,wages)*...
            (1/(1-param.beta*(1-param.delta)));

        iter = 1;    

        while err>tol

            Q_old = Q; 
            % 3. Compute the value of being single for each education level
            %    and sex
            [VS_f, VS_m] = val_single(param,expect,wages,Q);

            % 4. For each pair of female-male education levels, compute who 
            %    has the higher value of being single and record this value
            VS = max(ones(szP_f,szP_m).*VS_f,ones(szP_f,szP_m).*VS_m');
            
            % 5. Compute a new Q by equalizing the value of being married 
            %    to the value of being single for each pair of female-male 
            %    education levels. Repeat 3-5 until convergence.
            Q = (VS - VM)*(1-param.beta*(1-param.delta));
            err_report = sum(sum(abs(Q-Q_old)));

            fprintf('Iteration %d done, error was %g\n', iter, err_report);
            
            err = err_report;
            iter = iter+1;
            
        end
        
        [VS_f, VS_m] = val_single(param,expect,wages,Q);
        
        eps = exp(-10);
        ID = (ones(szP_f,szP_m).*VS_f>(ones(szP_f,szP_m).*VS_m'+eps))-...
             (ones(szP_f,szP_m).*VS_m'>(ones(szP_f,szP_m).*VS_f+eps));
               
        MS = (expect.P_f*expect.P_m').*(1-normcdf(Q,param.MU))/...
              sum(sum((expect.P_f*expect.P_m').*(1-normcdf(Q,param.MU))));
        
    else
        
        error_text_1 = ['Expected distribution of types and wage '...
            'vectors must have the same dimensions for each sex'];
        error(error_text_1);
        
    end
    
else
    
    error_text_2 = ['Expected distribution of types and wages must '...
        'be vectors'];
    error(error_text_2);
    
end



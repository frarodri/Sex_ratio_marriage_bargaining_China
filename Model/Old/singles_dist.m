%{
Computes the sex ratio and distribution of types among singles given a 
a matrix of cut-off qualities for marriage that is derived from the policy
function of agents in the PREE  
%}

function [theta_S, P_f, P_m] = singles_dist(Q,param)

if isscalar(param.theta)
    
    if size(Q,1)==length(param.P_f) && size(Q,2)==length(param.P_m)
        
        % Initialize algorithm by setting sex ratio to 1 and distributions
        % to uniform
        
        theta_S = 1;
        P_f = ones(length(param.P_f),1)/length(param.P_f);
        P_m = ones(length(param.P_m),1)/length(param.P_m);
        
        err = 1; % Inititalize error
        tol = 1*exp(-10); % Tolerance
        
        while err>tol
            
            theta_S_old = theta_S;
            P_f_old = P_f;
            P_m_old = P_m;
            
            [pi_f,pi_m] = match_prob(param,theta_S);
            
            % Calculate marriage rates
            MR_f = pi_f*(1-normcdf(Q,param.MU))*P_m; 
            MR_m = pi_m*(1-normcdf(Q,param.MU))'*P_f; 
            
            % Compute measure of single agents of each type
            M_S_f = param.P_f./(param.delta+MR_f-param.delta*MR_f); 
            M_S_m = param.P_m*param.theta./...
                    (param.delta+MR_m-param.delta*MR_m); 
            
            % Update sex ratio and distributions
            theta_S = sum(M_S_m)/sum(M_S_f);
            P_f = M_S_f/sum(M_S_f);
            P_m = M_S_m/sum(M_S_m);
            
            % Compute distance between previous and current iteration
            err = abs(theta_S-theta_S_old)+sum(abs(P_f-P_f_old))+...
                  sum(abs(P_m-P_m_old));
            
        end
                
    else
        
        error_text_1 = ['Number of rows and columns in cut-off quality '...
            'matrix must be equal to length of probability '...
            'distribution over types for females and males, respectively'];
        error(error_text_1);
        
    end
    
else
    
    error_text_2 = ['Sex ratio in the general population must be a '...
            'positive scalar'];
        error(error_text_2);
    
end


        
        
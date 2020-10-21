%{
Computes the steady state sex ratio and distribution of types among singles 
that would eventually arise given the matrix of marriage reservation match 
qualities Q that is derived from the policy function of agents in the PREE  
%}

function [theta_S, SPf, SPm] = dist_singles(param,Q)

if isscalar(param.theta) && param.theta>0
    
    if size(Q,1)==length(param.Pf) && size(Q,2)==length(param.Pm)
        
        % Initialize algorithm by setting sex ratio to 1 and distributions
        % to uniform
        
        theta_S = 1;
        SPf = ones(length(param.Pf),1)/length(param.Pf);
        SPm = ones(length(param.Pm),1)/length(param.Pm);
        
        % Calculate marriage probabilities among different couples types
        Mprob = (1-normcdf(Q,param.MU));
        
        current_error = 1; % Inititalize error
        tol = 1*exp(-10); % Tolerance
        
        while current_error>tol
            
            theta_S_prior = theta_S;
            SPf_prior = SPf;
            SPm_prior = SPm;
            
            [pi_f,pi_m] = match_prob(param,theta_S);
            
            % Calculate marriage rates
            MRf = pi_f*Mprob*SPm; 
            MRm = pi_m*Mprob'*SPf;
            
            % Compute measure of single agents of each type
            Sf = param.Pf./(param.delta+MRf-param.delta*MRf); 
            Sm = param.Pm*param.theta./...
                    (param.delta+MRm-param.delta*MRm);
                
            % Update sex ratio and distributions
            theta_S = sum(Sm)/sum(Sf);
            SPf = Sf/sum(Sf);
            SPm = Sm/sum(Sm);
            
            % Compute distance between previous and current iteration
            current_error = abs(theta_S-theta_S_prior)+...
                            sum(abs(SPf-SPf_prior))+...
                            sum(abs(SPm-SPm_prior));
            
        end
        
    else
        
        error_text_1 = ['Number of rows and columns in reservation '...
            'quality matrix must be equal to length of probability '...
            'distribution over types for females and males, respectively'];
        error(error_text_1);
        
    end
    
else
    
    error_text_2 = ['Sex ratio among new entrants must be a '...
            'positive scalar'];
    error(error_text_2);
    
end
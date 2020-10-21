%{
Calculates the matching probabilities for females and males given 
matching efficiency and elasticity parameters (param.A and param.alpha), 
and the expected sex ratio among singles theta_S
%}

function [pi_f , pi_m] = match_prob(param,theta_S)

if isscalar(param.A) && param.A>=0
    
    if isscalar(param.alpha) && param.alpha>=0 && param.alpha<=1
        
        if isscalar(theta_S) && theta_S>0
            
            X_m = [param.A*(1/theta_S)^(1-param.alpha) 1 1/theta_S];
            X_f = [param.A*theta_S^param.alpha theta_S 1];
            pi_f = min(X_f);
            pi_m = min(X_m);
            
        else
            
            error_text_1 = ['Male to female ratio among singles must '...
                'be a positive scalar'];
            error(error_text_1)
            
        end
        
    else
        
        error_text_2 = ['Matching function elasticity must be a scalar '...
            'between zero and one'];
        error(error_text_2)
        
    end
    
else
    
    error_text_3 = ['Matching function efficiency parameter must be a '...
        'non-negative scalar'];
    error(error_text_3)
    
 end


   
    


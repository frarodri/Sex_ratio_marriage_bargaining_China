function [theta_S,Q,PWf] = sseq(param,wages)
                 
current_err = 1;
tol = 10^(-5);
it = 1;

theta_S_prior=param.Pm*param.theta_0./param.Pf;

while current_err>tol
    
    [theta_S,Q,PWf] = theta_update(param,wages,theta_S_prior)
    current_err = sum(abs(theta_S-theta_S_prior))
    theta_S_prior = theta_S;
    
end
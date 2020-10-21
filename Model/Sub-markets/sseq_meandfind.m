function [theta_S,Q,PWf,MU] = ...
                      sseq_meandfind(param,wages,targ)

current_error = 1;
tol = 10^(-5);
it = 1;

theta_S_prior=param.theta_0*ones(size(param.Pf));

while current_error>tol
    
    [theta_S,Q,PWf,MU] = ...
        theta_update_meandfind(param,wages,targ,theta_S_prior)
    current_error = sum(abs(theta_S-theta_S_prior))
    theta_S_prior = theta_S;
    
end
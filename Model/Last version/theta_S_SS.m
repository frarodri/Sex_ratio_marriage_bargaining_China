%{
Computes the steady state sex ratio among singles in each submarket given a 
sex ratio among entrants, a reservation match quality and parameters
%}

function theta_S_new = theta_S_SS(param,theta_0,Q,theta_S_init)

theta_S_prior = theta_S_init;

current_err = 1;
tol = 10^(-5);

while current_err>tol
    
    [PIf, PIm] = match_prob(param,theta_S_prior);
    
    MRf = PIf.*(1-normcdf(diag(Q),diag(param.MU)));
    MRm = PIm.*(1-normcdf(diag(Q),diag(param.MU)));
    
    A = (1-(1-param.delta)*(1-MRf))./(1-(1-param.delta)*(1-MRm));
    
    theta_S_new = theta_0.*A;
    
    current_err = sum(abs(theta_S_prior-theta_S_new));
    
    theta_S_prior = theta_S_new;
    
end


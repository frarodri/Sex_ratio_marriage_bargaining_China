%{
Given a vector of different wages, computes the analytical solution to the 
problem faced each period by single agents
%}

function [C,L,G,H,EQ,N] = per_sol_singles(param,wage_vec)

if iscolumn(wage_vec) && all(wage_vec>0)
    
    XE = (1-param.alpha_g).*wage_vec/(param.alpha_g*param.pe);
    XG = param.A_g.*XE.^(1-param.alpha_g);   
    D = wage_vec./(param.A_g*param.alpha_g*XE.^(1-param.alpha_g));

    NUM_LM = (param.sigma_c^(1/param.sigma))*ones(size(wage_vec))+...
             wage_vec.*(param.sigma_l./wage_vec).^(1/param.sigma)+...
             (wage_vec.*(param.sigma_g./D).^(1/param.sigma))./XG+...
             param.pe*(XE./XG).*(param.sigma_g./D).^(1/param.sigma);

    LM = (NUM_LM./wage_vec).^param.sigma;

    C = (param.sigma_c./LM).^(1/param.sigma);
    L = (param.sigma_l./(LM.*wage_vec)).^(1/param.sigma);
    G = (param.sigma_g./(LM.*D)).^(1/param.sigma);
    H = G./XG;
    EQ = XE.*G./XG;
    N = 1-L-H;
    
else
    
    message = ['Vector of wages must be a column with strictly positive '...
        'entries'];
    error(message)
    
end
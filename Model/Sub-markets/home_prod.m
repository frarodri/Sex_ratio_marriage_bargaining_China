function G = home_prod(param,H,EQ)

if all(all(H>=0)) && all(all(EQ>=0))
    
    G = param.A_g*(EQ.^(1-param.alpha_g)).*H^param.alpha_g;
    
else
    
    message = ['Home production time and equipment matrices must be '...
        'non-negative'];
    error(message);

end


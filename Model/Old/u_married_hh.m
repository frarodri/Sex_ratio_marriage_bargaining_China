function u = u_married_hh(param,l_f)

    u = param.lambda*log(param.w_m+param.w_f*l_f-param.gamma*param.cbar)...
        +(1-param.lambda)*log(1-l_f);
    
end

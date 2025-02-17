function u = u_married(param,cf,cm,lf,lm,g,pwf)

pwm = 1-pwf;

if param.sigma == 1
    
    u = pwf*(param.sigma_c*log(cf)+param.sigma_l*log(lf))+...
        pwm*(param.sigma_c*log(cm)+param.sigma_l*log(lm))+...
        param.sigma_g*log(g);
    
else
    
    u = pwf*((param.sigma_c/(1-param.sigma))*cf^(1-param.sigma)+...
        (param.sigma_l/(1-param.sigma))*lf^(1-param.sigma))+...
        pwm*((param.sigma_c/(1-param.sigma))*cm^(1-param.sigma)+...
        (param.sigma_l/(1-param.sigma))*lm^(1-param.sigma))+...
        (param.sigma_g/(1-param.sigma))*g^(1-param.sigma);
    
end

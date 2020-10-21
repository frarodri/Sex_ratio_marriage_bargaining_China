%{
Computes the integrand that will be used to compute the expeted value of
marriage with Nash Bargaing for males
%}

function Im = integrand_nash_m(param,wages,VSf,VSm,mu,q)

pwf = nashbarg_pw(param,wages,VSf,VSm,q);

[~,UMm] = uflow_married(param,wages.f,wages.m,pwf,'quietly','true');

Mm = (UMm+q)/(1-param.beta*(1-param.delta));

Im = Mm*normpdf(q,mu);
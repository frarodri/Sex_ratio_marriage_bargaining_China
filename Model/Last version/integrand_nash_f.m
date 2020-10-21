%{
Computes the integrand that will be used to compute the expected value of
marriage with Nash Bargaing for females
%}

function If = integrand_nash_f(param,wages,VSf,VSm,mu,q)

pwf = nashbarg_pw(param,wages,VSf,VSm,q);

[UMf,~] = uflow_married(param,wages.f,wages.m,pwf,'quietly','true');

Mf = (UMf+q)/(1-param.beta*(1-param.delta));

If = Mf*normpdf(q,mu);


%%
% Test submkt_pree

%%
% Test parameters

param.beta = 0.96; 
param.sigma = 1.5; 
param.alpha_x = 1/2;
param.delta = 0.0203;
param.pe = 1;
param.A_g = 1;
param.alpha_g = 0.95;
param.eta = 0.33;
param.eta_f = 0.5470;
param.A_x = 1;
param.sigma_c = 0.3430;
param.sigma_l = 0.6180;
param.sigma_g = 0.039;
param.singcostf = 1.3;

param.Pf = ones(5,1)/5;
param.Pm = ones(5,1)/5;

wage_m = 1;
gender_wage_ratio = 0.75;
wage_f = gender_wage_ratio*wage_m;

theta_S = 10;

mu = 0;

%%
% Call submkt_pree

[qr,pw,VSf,VSm] = ...
    submkt_pree(param,theta_S,wage_f,wage_m,mu,'quietly','true');

[cf,cm,lf,lm,g,hf,hm,h,eq,nf,nm] = ...
                      per_sol_married(param,wages.f,wages.m,pw);

nm_hours = nm*118;
nf_hours = nf*118;

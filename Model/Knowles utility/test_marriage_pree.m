clear; clc;

param.sigma = 1.5;
param.sigma_c = 0.8;
param.sigma_l = 0.1;
param.sigma_g = 0.1;

param.eta = 0.33;
param.eta_f = 0.5;

param.A_g = 0.8;
param.alpha_g = 0.8;

param.pe = 1.5;

param.beta = 0.96;
param.delta = 1/45;

types = 5;
max_type = 2;

param.A_x = 0.8;
param.alpha_x = 0.5;
param.MU = zeros(types,types);

wages.m = linspace(1,max_type,types)';
wages.phi = 0.75;
wages.f = wages.phi*wages.m;

expect.theta_S = 1.8;
expect.Pf = ones(types,1)/types;
expect.Pm = expect.Pf;

PWf = 0.5*ones(types,types);

[Q,VSf,VSm,ID] = marriage_pree(param,expect,wages,PWf);

param.theta = 1.2;
param.Pf = ones(types,1)/types;
param.Pm = param.Pf;

init_expect.theta_S = param.theta;
init_expect.Pf = param.Pf;
init_expect.Pm = param.Pm;

[Q,theta_S,SPf,SPm] = ...
    marriage_sseq(param,wages,init_expect,PWf);

clear; clc;
types = 5;

param.sigma = 1.5;
param.sigma_c = 0.3;
param.sigma_l = 0.6;
param.sigma_g = 0.1;

param.eta = 0.3072;
param.eta_f = 0.6873;

param.A_g = 1;
param.alpha_g = 0.95;

param.pe = 1;

param.beta = 0.96;
param.delta = 1/49.230;

param.A_x = 0.8;
param.alpha_x = 0.5;
param.MU = zeros(types,types);

param.theta_zero = 1.0719;

param.Pf = [0.2553;0.2134;0.3399;0.1545;0.0369];
param.Pm = [0.1089;0.2192;0.4310;0.1962;0.0447];

wages.m = ones(5,1);
wages.phi = 0.7284;
wages.f = wages.phi*wages.m;

expect.theta_S = ones(types,1)*1.5;




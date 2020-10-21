clear ; close all; clc

%% 
% Baseline parameters

param.sigma = 1;
param.psi = -1.15;

% Parameters of the home production function
param.A_g = 5;
param.alpha_g = 0.95;

% Parameter for the price of home equipment
param.pe = 1;

% Wages
param.gender_wage_ratio = 0.7284;
param.gen_wage = 5;

% Wages
wage_m = param.gen_wage;
wage_f = param.gen_wage*param.gender_wage_ratio;

% Utility flow for singles
%  Females
param.sigma_c = 0.3526;
param.sigma_l = 0.6229;
param.sigma_g = 1-param.sigma_c-param.sigma_l;
USf = uflow_singles(param,wage_f);
% Males
param.sigma_c = 0.3635;
param.sigma_l = 0.6342;
param.sigma_g = 1-param.sigma_c-param.sigma_l;
USm = uflow_singles(param,wage_m);

% Parameters of the matching function
param.A_x = 1;
param.alpha_x = 1/2;

% Parameters of the effective discount rate
param.beta = 0.96;
param.delta = 1/45;

% Weights of the utility function for married people
param.sigma_c = 0.325;
param.sigma_l = 0.625;
param.sigma_g = 1-param.sigma_c-param.sigma_l;

% Parameters of the housework aggregator
param.eta_f = 0.55;
param.eta = 0.33;

% Mean of match quality draws
mu = 0;

% Sex ratio among singles
theta_S = 4;

% Starting point for algorithm
start_pw = 1/2;


%% 
% Compute the equilibrium

[qr,pwf,VSf,VSm] = ...
    submkt_pree(param,theta_S,wage_f,wage_m,mu,USf,USm,start_pw,...
    'Display','all','bargaining','Nash');

[cf,cm,lf,lm,g,hf,hm,h,eq,nf,nm] = ...
                      per_sol_married(param,wage_f,wage_m,pwf)



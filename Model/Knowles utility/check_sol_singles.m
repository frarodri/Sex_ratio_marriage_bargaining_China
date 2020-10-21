%**************************************************************************
% Checking the analytical solution for the single's problem               *
% Francisco Javier Rodríguez Román                                        * 
% Department of Economics                                                 *
% Universidad Carlos III de Madrid                                        *
%**************************************************************************

% home: G:\Mi unidad\China project\Model\Knowles utility
% uni: 

clear ; close all; clc
cd('G:\Mi unidad\China project\Model\Knowles utility')

% Parameters
param.sigma = 2;
param.sigma_c = 0.5;
param.sigma_l = 0.5;
param.sigma_g = 0.5;
param.A_g = 1;
param.alpha_g = 0.8;
param.pe = 0.5;

wage = 1;

% Numerical solution
A = [1,1,0];
b = 1;
Aeq = [];
beq = [];
x0 = [0.4,0.4,0.1];

lb = [0,0,0];
ub = [1,1,Inf];

[num_sol,ind_u] = ...
fmincon(@(controls)-obj_singles(param,wage,controls),x0,A,b,Aeq,beq,lb,ub);

num.c = bc_singles(param,wage,num_sol(1),num_sol(2),num_sol(3));
num.l = num_sol(2);
num.g = home_prod(param,num_sol(1),num_sol(3));
num.h = num_sol(1);
num.eq = num_sol(3);
num.n = 1-num.h-num.l;
num.ind_u = -ind_u; 

% Analytic solution

xe = (1-param.alpha_g)*wage/(param.alpha_g*param.pe);
xg = param.A_g*xe^(1-param.alpha_g);

D = wage/(param.A_g*param.alpha_g*xe^(1-param.alpha_g));

num_lagrange_m = param.sigma_c^(1/param.sigma)+...
                 wage*(param.sigma_l/wage)^(1/param.sigma)+...
                 (wage*(param.sigma_g/D)^(1/param.sigma))/xg+...
                 param.pe*(xe/xg)*(param.sigma_g/D)^(1/param.sigma);
             
lagrange_m = (num_lagrange_m/wage)^param.sigma;

an.c = (param.sigma_c/lagrange_m)^(1/param.sigma);
an.l = (param.sigma_l/(lagrange_m*wage))^(1/param.sigma);
an.g = (param.sigma_g/(lagrange_m*D))^(1/param.sigma);
an.h = an.g/xg;
an.eq = xe*an.g/xg;
an.n = 1-an.l-an.h;
an.ind_u = u_singles(param,an.c,an.l,an.g);                 


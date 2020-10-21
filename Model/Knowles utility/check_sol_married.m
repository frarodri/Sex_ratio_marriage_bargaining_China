%**************************************************************************
% Checking the analytical solution for the married couples' problem       *
% Francisco Javier Rodríguez Román                                        * 
% Department of Economics                                                 *
% Universidad Carlos III de Madrid                                        *
%**************************************************************************

% home: 
% uni: G:\Mi unidad\China project\Model\Knowles utility

clear ; close all; clc
cd('G:\Mi unidad\China project\Model\Knowles utility')

% Parameters
param.sigma = 2;
param.sigma_c = 0.5;
param.sigma_l = 0.5;
param.sigma_g = 0.5;
param.A_g = 0.8;
param.alpha_g = 0.8;
param.pe = 1.5;
param.eta = 0.33;
param.etaf = 0.52;

wages.f = 0.75;
wages.m = 2;
pwf = 0.5;
pwm = (1-pwf);

% Numerical solution
A = [0,0,1,0,1,0,0; 0,0,0,1,0,1,0;...
     1,1,wages.f,wages.m,wages.f,wages.m,param.pe];
b = [1;1;wages.f+wages.m];
Aeq = [];
beq = [];
x0 = [0.2,0.2,0.1,0.1,0.1,0.1,0.1];

lb = [0,0,0,0,0,0,0];
ub = [Inf,Inf,1,1,1,1,Inf];

[num_sol,ind_u] = ...
fmincon(@(ctrl)-obj_married(param,ctrl,pwf),x0,A,b,Aeq,beq,lb,ub);

num.cf = num_sol(1);    
num.cm = num_sol(2);
num.hf = num_sol(3);
num.hm = num_sol(4);                
num.lf = num_sol(5);
num.lm = num_sol(6);
num.h = h_married(param,num.hf,num.hm);
num.eq = num_sol(7);
num.g = home_prod(param,num.h,num.eq);
num.nf = 1-num.hf-num.lf;
num.nm = 1-num.hm-num.lm;
num.ind_u = -ind_u; 

% Analytical solution

xf = (param.etaf*wages.m/((1-param.etaf)*wages.f))^(1/param.eta);
xh = (param.etaf*xf^(1-param.eta)+1-param.etaf)^(1/(1-param.eta));
xe = ((1-param.alpha_g)*wages.m*xh^(1-param.eta))/...
    (param.alpha_g*param.pe*(1-param.etaf));
xg = param.A_g*(xe^(1-param.alpha_g))*xh^param.alpha_g;
D = (param.pe*(xe/xh)^param.alpha_g)/(param.A_g*(1-param.alpha_g));

num_lagrange_m = (pwf*param.sigma_c)^(1/param.sigma)+...
                 (pwm*param.sigma_c)^(1/param.sigma)+...
                 wages.f*(pwf*param.sigma_l/wages.f)^(1/param.sigma)+...
                 wages.m*(pwm*param.sigma_l/wages.m)^(1/param.sigma)+...
                 ((param.sigma_g/D)^(1/param.sigma))*...
                 (wages.f*xf/xg+wages.m/xg+param.pe*xe/xg);
         
lagrange_m = (num_lagrange_m/(wages.f+wages.m))^param.sigma;

an.cf = (pwf*param.sigma_c/lagrange_m)^(1/param.sigma);
an.cm = (pwm*param.sigma_c/lagrange_m)^(1/param.sigma);
g = (param.sigma_g/(lagrange_m*D))^(1/param.sigma);
an.hf = g*xf/xg;
an.hm = g/xg;
an.lf = (pwf*param.sigma_l/(lagrange_m*wages.f))^(1/param.sigma);
an.lm = (pwm*param.sigma_l/(lagrange_m*wages.m))^(1/param.sigma);
an.h = h_married(param,an.hf,an.hm);
an.eq = xe*g/xg;
an.g=g;
an.nf = 1-an.hf-an.lf;
an.nm = 1-an.hm-an.lm;
an.ind_u = u_married(param,an.cf,an.cm,an.lf,an.lm,an.g,pwf);











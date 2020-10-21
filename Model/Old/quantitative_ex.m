%**************************************************************************
% Quantitative exercise - second approach                                 *
% Francisco Javier Rodríguez Román                                        * 
% Department of Economics                                                 *
% Universidad Carlos III de Madrid                                        *
%**************************************************************************

% home: G:\My Drive\China project\Model
% uni: G:\Mi unidad\China project\Model

clear ; close all; clc
cd('G:\Mi unidad\China project\Model')

%%%%%%%%%%%%%%%
% Calibration %
%%%%%%%%%%%%%%%

% Parameters not calibrated

periods_year = 2;
param.beta = (0.96)^(1/periods_year); % Discount rate
life_expectancy = 45*periods_year; % Life expectancy at age 20
param.delta = 1/life_expectancy; % Death rate
param.alpha = 0.5; % Elasticity of matching function
param.gamma = 1;

% Data
% Sex ratio among singles
expect.theta_S = 628/470; 

% Distribution of education among singles
P_f = [70 ; 92 ; 203 ; 74 ; 31 ];
expect.P_f = P_f/sum(P_f);

P_m = [72 ; 137 ; 288 ; 104 ; 27 ];
expect.P_m = P_m/sum(P_m);

% Skill premium vectors
wages.gender_gap = .5599159/.7539827;

wages.OMEGA_m = ones(5,1);
wages.OMEGA_f = wages.gender_gap*wages.OMEGA_m;

% Female hours (relative to men)
data_mom.lf_bar = 2008.724/2248.348;

% Marital sorting
data_mom.delta = 1.633;
data_mom.h_f = 0.893;
data_mom.h_m = 0.618;

% Expected number of single periods
data_mom.exp_single_f = periods_year*3.13;
data_mom.exp_single_m = periods_year*5.17;

% Find the parameters that minimize the cost function

A = [];
b = []; 
Aeq = [];
beq = [];

x0 = [0.8,0.5,0,0.5];

lb = [0,0,-Inf,0.01];

ub = [1,0.7,Inf,1];

calibrated_params = ...
fmincon(@(param_to_calib)...
         cost_function(data_mom,param,expect,wages,param_to_calib),...
         x0,A,b,Aeq,beq,lb,ub);


[cost,model_mom] = ...
    cost_function(data_mom,param,expect,wages,calibrated_params);

param.lambda = calibrated_params(1);
param.cbar = calibrated_params(2);
param.kappa = calibrated_params(3);
param.A = calibrated_params(4);

param.MU = param.kappa*eye(length(expect.P_f));

% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Compute equilibrium for 2011 %
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Data

% Sex ratio among singles
expect.theta_S = 221/120; 

% Distribution of education among singles
P_f = [4 ; 3 ; 26 ; 11 ; 76 ];
expect.P_f = P_f/sum(P_f);

P_m = [3 ; 16 ; 65 ; 41 ; 96 ];
expect.P_m = P_m/sum(P_m);

% Skill premium vectors
real_wage_growth =  7.539237 /  2.234277;

OMEGA_m = [0; .2062788; .4133489; .6442774; 1.003563];
OMEGA_f = [0; .0828044 ; .2988064; .5541977; 1.004951];

OMEGA_m = 1+OMEGA_m;
OMEGA_f = 1+OMEGA_f;

gender_wage_gap_bottom = 1.965344 / 2.442528; 

wages.OMEGA_m = real_wage_growth*OMEGA_m;
wages.OMEGA_f = gender_wage_gap_bottom*real_wage_growth*OMEGA_f;

% Compute the equilibrium

[Q, VS_f, VS_m, ID, MS] = marriage_pree(param,expect,wages);

% Assortative matching
[delta_2011, h_f_2011, h_m_2011] = sort_cont_mat(MS);

% Labor supply
Lf_2011 = married_fls(param,wages);

lf_2011 = sum(sum(Lf_2011.*MS));

% Marriage rates/ periods spent single
[pi_f,pi_m] = match_prob(param,expect.theta_S);
                
mr_f_2011 = pi_f*expect.P_f'*(1-normcdf(Q,param.MU))*expect.P_m;
mr_m_2011 = pi_m*expect.P_m'*(1-normcdf(Q,param.MU))'*expect.P_f;
                
exp_single_f_2011 = 1/mr_f_2011;
exp_single_m_2011 = 1/mr_m_2011;

% Inequality
[ineq_inc_2011, ineq_cons_2011] = ...
    ineq(param,wages,Q,MS,expect.theta_S,expect.P_f,expect.P_m);

names = ...
    {'No schooling'; 'Primary'; 'Lower Middle'; 'Upper Middle'; 'College'};

Lf_2011_fig = figure('Color','White');
plot([1,2,3,4,5],Lf_2011')
legend(names,'Location','northeast')
legend('boxoff')
xlabel('Husband''s education')
ylabel('Labor supply (fraction of husband''s)')
set(gca,'xtick',[1:5],'xticklabel',names)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Exercise 1: sex ratio stays at 1991 level %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%5%%%%

% Sex ratio among singles
expect.theta_S = 628/470;

[Q, VS_f, VS_m, ID, MS] = marriage_pree(param,expect,wages);

% Assortative matching
[delta_ex1, h_f_ex1, h_m_ex1] = sort_cont_mat(MS);

% Labor supply
Lf_ex1 = married_fls(param,wages);

lf_ex1 = sum(sum(Lf_ex1.*MS));

% Marriage rates/ periods spent single
[pi_f,pi_m] = match_prob(param,expect.theta_S);
                
mr_f_ex1 = pi_f*expect.P_f'*(1-normcdf(Q,param.MU))*expect.P_m;
mr_m_ex1 = pi_m*expect.P_m'*(1-normcdf(Q,param.MU))'*expect.P_f;
                
exp_single_f_ex1 = 1/mr_f_ex1;
exp_single_m_ex1 = 1/mr_m_ex1;

% Inequality
[ineq_inc_ex1, ineq_cons_ex1] = ...
    ineq(param,wages,Q,MS,expect.theta_S,expect.P_f,expect.P_m);

% Check smoothness of objective function

% param.lambda = 0.5;
% param.cbar = 0.5;
% param.gamma = 1.5;
% param.kappa = 0.2;
% param.A = 0.5;
% 
% values = (0:0.01:0.74);
% obj_values = zeros(size(values));
% 
% for i = 1:length(values)
% 
%     param.cbar = values(i);
%     param.MU = param.kappa*eye(5); 
%     
%     obj_values(i) = obj_func(param,expect,wages,data_mom);
%     
% end
% 
% plot(values,obj_values)


% Looks smooth







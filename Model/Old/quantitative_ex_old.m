%**************************************************************************
% Quantitative exercise - first approach                                  *
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

periods_year = 4;
param.beta = (0.96)^(1/periods_year); % Discount rate
life_expectancy = 45*periods_year; % Life expectancy at age 20
param.delta = 1/life_expectancy; % Death rate
param.alpha = 0.5; % Elasticity of matching function

% Data
% Sex ratio among singles
param.theta_S = 628/470; 

% Distribution of education among singles
Z_f = [70 ; 92 ; 203 ; 74 ; 31 ];
Z_f = Z_f/sum(Z_f);

Z_m = [72 ; 137 ; 288 ; 104 ; 27 ];
Z_m = Z_m/sum(Z_m);

% Skill premium vectors
gender_wage_gap = .5599159/.7539827;

OMEGA_m = ones(5,1);
OMEGA_f =gender_wage_gap*OMEGA_m;

% Female hours (relative to men)
l_f_data = 2008.724/2248.348;

% Marital sorting
delta_data = 1.633;
h_f_data = 0.893;
h_m_data = 0.618;

% Expected number of single periods
exp_single_f_data = periods_year*3.13;
exp_single_m_data = periods_year*5.17;

% Pre-assigning matrices
min_lambda = (l_f_data+1/gender_wage_gap-1)/...
            (1+1/gender_wage_gap-1);
        
max_lambda = (l_f_data+1/gender_wage_gap)/...
            (1+1/gender_wage_gap);
eps = 0.001;
        
values_lambda = linspace(min_lambda+eps,max_lambda,101);

values_cbar = gender_wage_gap*(l_f_data./(1-values_lambda)...
    -values_lambda./(1-values_lambda)+1/gender_wage_gap);

% values_cbar = (l_f-values_lambda./(1-values_lambda)+1/gender_wage_gap)...
%     *gender_wage_gap;

values_A = linspace(0.1,1,37);
values_kappa = linspace(0,-0.5,21);

values_delta = ...
    zeros(length(values_lambda),length(values_A),length(values_kappa));
values_h_f = ...
    zeros(length(values_lambda),length(values_A),length(values_kappa));
values_h_m = ...
    zeros(length(values_lambda),length(values_A),length(values_kappa));

values_exp_single_f = ...
    zeros(length(values_lambda),length(values_A),length(values_kappa));
values_exp_single_m = ...
    zeros(length(values_lambda),length(values_A),length(values_kappa));

for i=1:length(values_lambda)
   
    for j=1:length(values_A)
       
        for k=1:length(values_kappa)
            
            if values_cbar(i)<0 || values_cbar(i)>gender_wage_gap
                
                values_delta(i,j,k) = NaN;
                values_h_f(i,j,k) = NaN;
                values_h_m(i,j,k) = NaN;
            
            else
               
                param.lambda = values_lambda(i);
                param.cbar = values_cbar(i);
                param.A = values_A(j);
                param.kappa = values_kappa(k);
                param.MU = param.kappa*(1-eye(length(Z_f)));
                
                [Q, VS_f, VS_m, ID, MS] = ...
                marriage_eq(param,OMEGA_f,OMEGA_m,Z_f,Z_m);
            
                [delta, h_f, h_m] = sort_cont_mat(MS);
                
                values_delta(i,j,k) = delta;
                values_h_f(i,j,k) = h_m;
                values_h_m(i,j,k) = h_f;
                
                [pi_f,pi_m] = match_prob(param);
                
                mr_f = pi_f*Z_f'*(1-normcdf(Q,param.MU))*Z_m;
                mr_m = pi_m*Z_f'*(1-normcdf(Q,param.MU))*Z_m;
                
                exp_single_f = 1/mr_f;
                exp_single_m = 1/mr_m;
                
                values_exp_single_f(i,j,k) = exp_single_f;
                values_exp_single_m(i,j,k) = exp_single_m;
                
            end
                
        end
        
    end
    
end

% values_error = (1/delta_data)*abs(values_delta-delta_data)+...
%                (1/h_f_data)*abs(values_h_f-h_f_data)+...
%                (1/h_m_data)*abs(values_h_m-h_m_data)+...
%                (1/exp_single_f_data)*...
%                abs(values_exp_single_f-exp_single_f_data)+...
%                (1/exp_single_m_data)*...
%                abs(values_exp_single_m-exp_single_m_data);

values_error = (1/delta_data)^2*abs(values_delta-delta_data)+...
               (1/exp_single_f_data)^2*...
               abs(values_exp_single_f-exp_single_f_data)+...
               (1/exp_single_m_data)^2*...
               abs(values_exp_single_m-exp_single_m_data);


[ error , I_kappa] = min(min(min(values_error)));
[ ~ , I_A] = min(min(values_error(:,:,I_kappa)));
[ ~ , I_lambda] = min(values_error(:,I_A,I_kappa));

param.lambda = values_lambda(I_lambda);
param.cbar = values_cbar(I_lambda);
param.A = values_A(I_A);
param.kappa = values_kappa(I_kappa);
param.MU = param.kappa*(1-eye(length(Z_f)));

delta_model_direct = values_delta(I_lambda,I_A,I_kappa);
h_f_model = values_h_f(I_lambda,I_A,I_kappa);
h_m_model = values_h_m(I_lambda,I_A,I_kappa);
exp_single_f_model_direct = values_exp_single_f(I_lambda,I_A,I_kappa);
exp_single_m_model_direct = values_exp_single_m(I_lambda,I_A,I_kappa);

[Q, VS_f, VS_m, ID, MS] = ...
                marriage_eq(param,OMEGA_f,OMEGA_m,Z_f,Z_m);
            
 % Assortative matching
[delta_model, ~, ~] = sort_cont_mat(MS);

% Labor supply
Lf_model = married_fls(param,OMEGA_f,OMEGA_m);

lf_model = sum(sum(Lf_model.*MS));

% Marriage rates/ periods spent single
[pi_f,pi_m] = match_prob(param);
                
mr_f_model = pi_f*Z_f'*(1-normcdf(Q,param.MU))*Z_m;
mr_m_model = pi_m*Z_f'*(1-normcdf(Q,param.MU))*Z_m;
                
exp_single_f_model = 1/mr_f_model;
exp_single_m_model = 1/mr_m_model;

% Inequality
[ineq_inc_model, ineq_cons_model] = ...
    ineq(param,Q,MS,Z_f,Z_m,OMEGA_f,OMEGA_m);
            
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Compute equilibrium for 2011 %
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Data

% Sex ratio among singles
param.theta_S = 221/120; 

% Distribution of education among singles
Z_f = [4 ; 3 ; 26 ; 11 ; 76 ];
Z_f = Z_f/sum(Z_f);

Z_m = [3 ; 16 ; 65 ; 41 ; 96 ];
Z_m = Z_m/sum(Z_m);

% Skill premium vectors
real_wage_growth =  7.539237 /  2.234277;

OMEGA_m = [0; .2062788; .4133489; .6442774; 1.003563];
OMEGA_f = [0; .0828044 ; .2988064; .5541977; 1.004951];

OMEGA_m = 1+OMEGA_m;
OMEGA_f = 1+OMEGA_f;

gender_wage_gap_bottom = 1.965344 / 2.442528; 

OMEGA_m = real_wage_growth*OMEGA_m;
OMEGA_f = gender_wage_gap_bottom*real_wage_growth*OMEGA_f;

% Compute the equilibrium

[Q, VS_f, VS_m, ID, MS] = ...
    marriage_eq(param,OMEGA_f,OMEGA_m,Z_f,Z_m);

% Assortative matching
[delta_2011, h_f_2011, h_m_2011] = sort_cont_mat(MS);

% Labor supply
Lf_2011 = married_fls(param,OMEGA_f,OMEGA_m);

lf_2011 = sum(sum(Lf_2011.*MS));

% Marriage rates/ periods spent single
[pi_f,pi_m] = match_prob(param);
                
mr_f_2011 = pi_f*Z_f'*(1-normcdf(Q,param.MU))*Z_m;
mr_m_2011 = pi_m*Z_f'*(1-normcdf(Q,param.MU))*Z_m;
                
exp_single_f_2011 = 1/mr_f_2011;
exp_single_m_2011 = 1/mr_m_2011;

% Inequality
[ineq_inc_2011, ineq_cons_2011] = ineq(param,Q,MS,Z_f,Z_m,OMEGA_f,OMEGA_m);

names = {'No schooling'; 'Primary'; 'Lower Middle'; 'Upper Middle'; 'College'};

Lf_2011_fig = figure('Color','White');
plot([1,2,3,4,5],Lf_2011')
legend(names,'Location','northeast')
legend('boxoff')
xlabel('Husband''s education')
ylabel('Labor supply (fraction of husband''s)')
set(gca,'xtick',[1:5],'xticklabel',names)

% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% % Exercise 1: sex ratio stays at 1991 level %
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%5%%%%
% 
% % Sex ratio among singles
% param.theta_S = 628/470;
% 
% [Q, VS_f, VS_m, ID, MS] = ...
%     marriage_eq(param,OMEGA_f,OMEGA_m,Z_f,Z_m);
% 
% % Assortative matching
% [delta_ex1, h_f_ex1, h_m_ex1] = sort_cont_mat(MS);
% 
% % Labor supply
% Lf_ex1 = married_fls(param,OMEGA_f,OMEGA_m);
% 
% lf_ex1 = sum(sum(Lf_ex1.*MS));
% 
% % Marriage rates/ periods spent single
% [pi_f,pi_m] = match_prob(param);
%                 
% mr_f_ex1 = pi_f*Z_f'*(1-normcdf(Q,param.MU))*Z_m;
% mr_m_ex1 = pi_m*Z_f'*(1-normcdf(Q,param.MU))*Z_m;
%                 
% exp_single_f_ex1 = 1/mr_f_ex1;
% exp_single_m_ex1 = 1/mr_m_ex1;
% 
% % Inequality
% [ineq_inc_ex1, ineq_cons_ex1] = ineq(param,Q,MS,Z_f,Z_m,OMEGA_f,OMEGA_m);
% 
% 
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% % Exercise 2: sex ratio increases further to 2.5 %
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 
% % Sex ratio among singles
% param.theta_S = 2.5;
% 
% [Q, VS_f, VS_m, ID, MS] = ...
%     marriage_eq(param,OMEGA_f,OMEGA_m,Z_f,Z_m);
% 
% % Assortative matching
% [delta_ex2, h_f_ex2, h_m_ex2] = sort_cont_mat(MS);
% 
% % Labor supply
% Lf_ex2 = married_fls(param,OMEGA_f,OMEGA_m);
% 
% lf_ex2 = sum(sum(Lf_ex2.*MS));
% 
% % Marriage rates/ periods spent single
% [pi_f,pi_m] = match_prob(param);
%                 
% mr_f_ex2 = pi_f*Z_f'*(1-normcdf(Q,param.MU))*Z_m;
% mr_m_ex2 = pi_m*Z_f'*(1-normcdf(Q,param.MU))*Z_m;
%                 
% exp_single_f_ex2 = 1/mr_f_ex2;
% exp_single_m_ex2 = 1/mr_m_ex2;
% 
% % Inequality
% [ineq_inc_ex2, ineq_cons_ex2] = ineq(param,Q,MS,Z_f,Z_m,OMEGA_f,OMEGA_m);
% 
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% % Exercise 3: no gender wage gap %
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 
% % Sex ratio among singles
% param.theta_S = 221/120;
% 
% % Wages
% OMEGA_f = OMEGA_m;
% 
% [Q, VS_f, VS_m, ID, MS] = ...
%     marriage_eq(param,OMEGA_f,OMEGA_m,Z_f,Z_m);
% 
% % Assortative matching
% [delta_ex3, h_f_ex3, h_m_ex3] = sort_cont_mat(MS);
% 
% % Labor supply
% Lf_ex3 = married_fls(param,OMEGA_f,OMEGA_m);
% 
% lf_ex3 = sum(sum(Lf_ex3.*MS));
% 
% % Marriage rates/ periods spent single
% [pi_f,pi_m] = match_prob(param);
%                 
% mr_f_ex3 = pi_f*Z_f'*(1-normcdf(Q,param.MU))*Z_m;
% mr_m_ex3 = pi_m*Z_f'*(1-normcdf(Q,param.MU))*Z_m;
%                 
% exp_single_f_ex3 = 1/mr_f_ex3;
% exp_single_m_ex3 = 1/mr_m_ex3;
% 
% % Inequality
% [ineq_inc_ex3, ineq_cons_ex3] = ineq(param,Q,MS,Z_f,Z_m,OMEGA_f,OMEGA_m);


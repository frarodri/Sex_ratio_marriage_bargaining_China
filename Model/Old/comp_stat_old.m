%**************************************************************************
% Comparative statics                                                     *
% Francisco Javier Rodríguez Román                                        * 
% Department of Economics                                                 *
% Universidad Carlos III de Madrid                                        *
%**************************************************************************

% home: G:\My Drive\China project\Model
% uni: G:\Mi unidad\China project\Model

clear ; close all; clc
cd('G:\Mi unidad\China project\Model')
changing_qty = 'theta_S';
changing_param = 'lambda';

% Baseline paramaters

periods_year = 1;
param.beta = (0.96)^(1/periods_year); % Discount rate
life_expectancy = 45*periods_year; % Life expectancy at age 20
param.delta = 1/life_expectancy; % Death rate
param.A = 0.5; % Efficiency of matching function
param.alpha = 0.5; % Elasticity of matching function
param.MU = 0*eye(5); % Mean of the distribution of match quality by 
                       % partners' education 
param.cbar = 0; % Household upkeep cost
param.lambda = 1; % Weight of market goods consumption in household 
                  % utility
param.theta_S = 1.34; % Sex ratio among singles
top = 2;
param.skill_premium = linspace(1,top,5)'; % Skill premium
param.gender_wage_gap = 1; % Gender wage gap

% From the data

% Distribution of education among singles
% Z_f = [0.0932 ; 0.087 ; 0.3478 ; 0.2981 ; 0.1739];
% Z_m = [0.1019 ; 0.088 ; 0.4167 ; 0.2730 ; 0.1204];

Z_f = ones(5,1)/5;
Z_m = Z_f;

% Skill premium vectors
OMEGA_m = param.skill_premium;
OMEGA_f = param.gender_wage_gap*OMEGA_m;

% Changing the sex ratio among singles
values_theta_S = linspace(0.4,2.5,201);
VS_f_theta_S = zeros(size(Z_f,1),size(values_theta_S,2)); 
VS_m_theta_S = zeros(size(Z_m,1),size(values_theta_S,2));
id_theta_S = zeros(size(values_theta_S));
delta_theta_S = zeros(size(values_theta_S));
hyper_f_theta_S = zeros(size(values_theta_S));
hyper_m_theta_S = zeros(size(values_theta_S));
delta_theta_S2 = zeros(size(values_theta_S));
hyper_f_theta_S2 = zeros(size(values_theta_S));
hyper_m_theta_S2 = zeros(size(values_theta_S));
MR_f_theta_S = zeros(size(Z_f,1),size(values_theta_S,2)); 
MR_m_theta_S = zeros(size(Z_m,1),size(values_theta_S,2));
mr_f_theta_S = zeros(size(values_theta_S));
mr_m_theta_S = zeros(size(values_theta_S));
exp_single_f_theta_S = zeros(size(values_theta_S));
exp_single_m_theta_S = zeros(size(values_theta_S));
lf_bar_theta_S = zeros(size(values_theta_S));
pi_f_theta_S = zeros(size(values_theta_S));
pi_m_theta_S = zeros(size(values_theta_S));

for i=1:length(values_theta_S)
    
    param.theta_S = values_theta_S(i);
    [Q, VS_f, VS_m, ID, MS] = ...
    marriage_eq(param,OMEGA_f,OMEGA_m,Z_f,Z_m);

    % Value of being single
    VS_f_theta_S(:,i) = VS_f;
    VS_m_theta_S(:,i) = VS_m;
    
    % Percentage of matches in which female is binding party
    id_theta_S(i) = sum(sum(Z_f*Z_m'.*(ID==1)));
    
    % Marital sorting
    [delta, h_f, h_m] = sort_cont_mat(MS);
    delta_theta_S(i) = delta;
    hyper_f_theta_S(i) = h_f;
    hyper_m_theta_S(i) = h_m;
    
    [delta2, h_f2, h_m2] = sort_cont_mat_2(MS,Z_f,Z_m);
    delta_theta_S2(i) = delta2;
    hyper_f_theta_S2(i) = h_f2;
    hyper_m_theta_S2(i) = h_m2;
    
    % Marriage rates
    [pi_f,pi_m] = match_prob(param);
    MR_f_theta_S(:,i) = ...
        sum(pi_f*Z_f*Z_m'.*(1-normcdf(Q,param.MU)),2)./Z_f;
    MR_m_theta_S(:,i) = ...
        sum(pi_m*Z_f*Z_m'.*(1-normcdf(Q,param.MU)),1)'./Z_m;


    mr_f_theta_S(i) = pi_f*Z_f'*(1-normcdf(Q,param.MU))*Z_m;
    mr_m_theta_S(i) = pi_m*Z_f'*(1-normcdf(Q,param.MU))*Z_m;
    
    % Expected number of periods spent in single's pool
    exp_single_f_theta_S(i) = 1/mr_f_theta_S(i);
    exp_single_m_theta_S(i) = 1/mr_m_theta_S(i);
    
    % Aggregate labor supply of married women relative to men
    Lf = married_fls(param,OMEGA_f,OMEGA_m);
    lf_bar_theta_S(i) = sum(sum(Lf.*MS));
    
    % Match probabilities
    [pi_f_theta_S(i),pi_m_theta_S(i)] = match_prob(param);
    
end

VS_f_theta_S_fig = figure('Color','White');
plot(values_theta_S,VS_f_theta_S);
legend({'No schooling','Primary','Lower middle','Upper middle','College'}...
,'Location','northeast');
legend('boxoff')
xlabel('\theta_S')
ylabel('Value of being single')

VS_m_theta_S_fig = figure('Color','White');
plot(values_theta_S,VS_m_theta_S);
legend({'No schooling','Primary','Lower middle','Upper middle','College'}...
,'Location','northeast');
legend('boxoff')
xlabel('\theta_S')
ylabel('Value of being single')

sorting_theta_S_fig = figure('Color','White');
plot(values_theta_S,delta_theta_S);
legend('boxoff')
xlabel('\theta_S')
ylabel('\delta')

hypergamy_theta_S_fig = figure('Color','White');
plot(values_theta_S,[hyper_f_theta_S; hyper_m_theta_S])
legend({'Female', 'Male'},'Location','northeast')
legend('boxoff')
xlabel('\theta_S')
ylabel('Hypergamy measure')

sorting_theta_S2_fig = figure('Color','White');
plot(values_theta_S,delta_theta_S2);
legend('boxoff')
xlabel('\theta_S')
ylabel('\delta')

hypergamy_theta_S2_fig = figure('Color','White');
plot(values_theta_S,[hyper_f_theta_S2; hyper_m_theta_S2])
legend({'Female', 'Male'},'Location','northeast')
legend('boxoff')
xlabel('\theta_S')
ylabel('Hypergamy measure')

id_theta_S_fig = figure('Color','White');
plot(values_theta_S,id_theta_S);
xlabel('\theta_S')
ylabel('Fraction of matches female is pivotal')

MR_f_theta_S_fig = figure('Color','White');
plot(values_theta_S,MR_f_theta_S);
legend({'No schooling','Primary','Lower middle','Upper middle','College'}...
,'Location','northwest');
legend('boxoff')
xlabel('\theta_S')
ylabel('Marriage rate per period')

MR_m_theta_S_fig = figure('Color','White');
plot(values_theta_S,MR_m_theta_S);
legend({'No schooling','Primary','Lower middle','Upper middle','College'}...
,'Location','northeast');
legend('boxoff')
xlabel('\theta_S')
ylabel('Marriage rate per period')

marr_rates_theta_S_fig = figure('Color','White');
plot(values_theta_S,[mr_f_theta_S; mr_m_theta_S]);
legend({'Females','Males'},'Location','northwest')
legend('boxoff')
xlabel('\theta_S')
ylabel('Marriage rate per period')

match_prob_theta_S_fig = figure('Color','White');
plot(values_theta_S,[pi_f_theta_S; pi_m_theta_S]);
legend({'Females','Males'},'Location','northwest')
legend('boxoff')
xlabel('\theta_S')
ylabel('Matching probability')

exp_single_A_fig = figure('Color','White');
plot(values_theta_S,[exp_single_f_theta_S; exp_single_m_theta_S]);
legend({'Females','Males'},'Location','northwest')
legend('boxoff')
xlabel('\theta_S')
ylabel('Expected number of periods spent single')

lf_theta_S_fig = figure('Color','White');
plot(values_theta_S,lf_bar_theta_S);
legend('boxoff')
xlabel('\theta_S')
ylabel('Married female market hours as fraction of husbands')

VS_theta_S_fig = figure('Color','White');
plot(values_theta_S,[VS_f_theta_S;VS_m_theta_S]);
xlabel('\theta_S')
ylabel('Value of being single')
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

%{
Scenario 1: no inequality of any kind (no skill premium or gender gap), in
which case the distribution of types shouldn't matter 
%}

% Parameters
types_f = 5; % Number of female types
types_m = 5; % Number of male types
periods_year = 1; % How many times during the year are we going to allow 
                  % matches to occur
param.beta = (0.96)^(1/periods_year); % Discount rate
life_expectancy = 45*periods_year; % Life expectancy in the first period
param.delta = 1/life_expectancy; % Death rate
param.A = 1; % Efficiency of matching
param.alpha = 0.5; % Elasticity of matching function
param.cbar = 0; % Household upkeep cost
param.lambda = 1; % Weight of market goods consumption in household's 
                  % utility
param.theta_S = 1; % Sex ratio among singles
param.MU = zeros(types_f,types_m); % Mean of match quality for each 
                                   % female-male type combination
param.gender_wage_gap = 1; % Gender wage gap
                                   
% Distributions of income and types
OMEGA_m = ones(types_m,1);
OMEGA_f = param.gender_wage_gap*OMEGA_m;

Z_f = ones(types_f,1)/types_f;
Z_m = ones(types_m,1)/types_m;

% Changing efficiency of matching
vals_A = linspace(0.01,1,100);
id_A = zeros(size(vals_A));
delta_A = zeros(size(vals_A));
hyper_f_A = zeros(size(vals_A));
hyper_m_A = zeros(size(vals_A));
MR_f_A = zeros(size(Z_f,1),size(vals_A,2)); 
MR_m_A = zeros(size(Z_m,1),size(vals_A,2));
mr_f_A = zeros(size(vals_A));
mr_m_A = zeros(size(vals_A));
exp_single_f_A = zeros(size(vals_A));
exp_single_m_A = zeros(size(vals_A));
lf_bar_A = zeros(size(vals_A));
pi_f_A = zeros(size(vals_A));
pi_m_A = zeros(size(vals_A));

for i=1:length(vals_A)
    
    param.A = vals_A(i);
    [Q, ID] = marriage_eq(param,OMEGA_f,OMEGA_m,Z_f,Z_m);
    
    % Percentage of matches in which female is binding party
    id_A(i) = sum(sum(Z_f*Z_m'.*(ID>=0)));
    
    % Marital sorting
    MS = (Z_f*Z_m'.*(1-normcdf(Q,param.MU)))...
        /sum(sum(Z_f*Z_m'.*(1-normcdf(Q,param.MU))));
    delta_A(i) = sort_cont_mat(MS);
    [hyper_f_A(i), hyper_m_A(i)] = hypergamy(MS);
    
    % Marriage rates
    [pi_f,pi_m] = match_prob(param);
    MR_f_A(:,i) = ...
        sum(pi_f*Z_f*Z_m'.*(1-normcdf(Q,param.MU)),2)./Z_f;
    MR_m_A(:,i) = ...
        sum(pi_m*Z_f*Z_m'.*(1-normcdf(Q,param.MU)),1)'./Z_m;


    mr_f_A(i) = pi_f*Z_f'*(1-normcdf(Q,param.MU))*Z_m;
    mr_m_A(i) = pi_m*Z_f'*(1-normcdf(Q,param.MU))*Z_m;
    % Expected number of periods spent in single's pool
    exp_single_f_A(i) = 1/mr_f_A(i);
    exp_single_m_A(i) = 1/mr_m_A(i);
    
    % Aggregate labor supply of married women relative to men
    Lf = married_fls(param,OMEGA_f,OMEGA_m);
    lf_bar_A(i) = sum(sum(Lf.*MS));
    
    % Match probabilities
    [pi_f_A(i),pi_m_A(i)] = match_prob(param);
    
end

% Graphs

sorting_A_fig = figure('Color','White');
plot(vals_A,delta_A);
legend('boxoff')
xlabel('Matching efficiency (A)')
yticks(1)
ylabel('\delta')

hypergamy_A_fig = figure('Color','White');
plot(vals_A,[hyper_f_A; hyper_m_A])
legend({'Female', 'Male'},'Location','northeast')
legend('boxoff')
xlabel('Matching efficiency (A)')
yticks(1)
ylabel('Hypergamy measure')

id_A_fig = figure('Color','White');
plot(vals_A,id_A);
xlabel('Matching efficiency (A)')
ylabel('Fraction of matches female is pivotal')

MR_f_A_fig = figure('Color','White');
plot(vals_A,MR_f_A);
legend({'No schooling','Primary','Lower middle','Upper middle','College'}...
,'Location','northwest');
legend('boxoff')
xlabel('Matching efficiency (A)')
ylabel('Marriage rate per period')

MR_m_A_fig = figure('Color','White');
plot(vals_A,MR_m_A);
legend({'No schooling','Primary','Lower middle','Upper middle','College'}...
,'Location','northwest');
legend('boxoff')
xlabel('Matching efficiency (A)')
ylabel('Marriage rate per period')

marr_rates_A_fig = figure('Color','White');
plot(vals_A,[mr_f_A; mr_m_A]);
legend({'Females','Males'},'Location','northwest')
legend('boxoff')
xlabel('Matching efficiency (A)')
ylabel('Marriage rate per period')

match_prob_A_fig = figure('Color','White');
plot(vals_A,[pi_f_A; pi_m_A]);
legend({'Females','Males'},'Location','northwest')
legend('boxoff')
xlabel('Matching efficiency (A)')
ylabel('Matching probability')

exp_single_A_fig = figure('Color','White');
plot(vals_A,[exp_single_f_A; exp_single_m_A]);
legend({'Females','Males'},'Location','northwest')
legend('boxoff')
xlabel('Matching efficiency (A)')
ylabel('Expected number of periods spent single')

lf_A_fig = figure('Color','White');
plot(vals_A,lf_bar_A);
legend('boxoff')
xlabel('Matching efficiency (A)')
yticks(1)
ylabel('Married female market hours as fraction of husbands')

% Changing the weight of market goods consumption in household's utility
param.A = 0.5;

vals_lambda = linspace(0,1,101);
id_lambda = zeros(size(vals_lambda));
delta_lambda = zeros(size(vals_lambda));
hyper_f_lambda = zeros(size(vals_lambda));
hyper_m_lambda = zeros(size(vals_lambda));
MR_f_lambda = zeros(size(Z_f,1),size(vals_lambda,2)); 
MR_m_lambda = zeros(size(Z_m,1),size(vals_lambda,2));
mr_f_lambda = zeros(size(vals_lambda));
mr_m_lambda = zeros(size(vals_lambda));
exp_single_f_lambda = zeros(size(vals_lambda));
exp_single_m_lambda = zeros(size(vals_lambda));
lf_bar_lambda = zeros(size(vals_lambda));
pi_f_lambda = zeros(size(vals_lambda));
pi_m_lambda = zeros(size(vals_lambda));

for i=1:length(vals_lambda)
    
    param.lambda = vals_lambda(i);
    [Q, ID] = marriage_eq(param,OMEGA_f,OMEGA_m,Z_f,Z_m);
    
    % Percentage of matches in which female is binding party
    id_lambda(i) = sum(sum(Z_f*Z_m'.*(ID>=0)));
    
    % Marital sorting
    MS = (Z_f*Z_m'.*(1-normcdf(Q,param.MU)))...
        /sum(sum(Z_f*Z_m'.*(1-normcdf(Q,param.MU))));
    delta_lambda(i) = sort_cont_mat(MS);
    [hyper_f_lambda(i), hyper_m_lambda(i)] = hypergamy(MS);
    
    % Marriage rates
    [pi_f,pi_m] = match_prob(param);
    MR_f_lambda(:,i) = ...
        sum(pi_f*Z_f*Z_m'.*(1-normcdf(Q,param.MU)),2)./Z_f;
    MR_m_lambda(:,i) = ...
        sum(pi_m*Z_f*Z_m'.*(1-normcdf(Q,param.MU)),1)'./Z_m;


    mr_f_lambda(i) = pi_f*Z_f'*(1-normcdf(Q,param.MU))*Z_m;
    mr_m_lambda(i) = pi_m*Z_f'*(1-normcdf(Q,param.MU))*Z_m;
    % Expected number of periods spent in single's pool
    exp_single_f_lambda(i) = 1/mr_f_lambda(i);
    exp_single_m_lambda(i) = 1/mr_m_lambda(i);
    
    % Aggregate labor supply of married women relative to men
    Lf = married_fls(param,OMEGA_f,OMEGA_m);
    lf_bar_lambda(i) = sum(sum(Lf.*MS));
    
    % Match probabilities
    [pi_f_lambda(i),pi_m_lambda(i)] = match_prob(param);
    
end

% Graphs

sorting_lambda_fig = figure('Color','White');
plot(vals_lambda,delta_lambda);
legend('boxoff')
xlabel('\lambda')
yticks(1)
ylabel('\delta')

hypergamy_lambda_fig = figure('Color','White');
plot(vals_lambda,[hyper_f_lambda; hyper_m_lambda])
legend({'Female', 'Male'},'Location','northeast')
legend('boxoff')
xlabel('\lambda')
yticks(1)
ylabel('Hypergamy measure')

id_lambda_fig = figure('Color','White');
plot(vals_lambda,id_lambda);
xlabel('\lambda')
ylabel('Fraction of matches female is pivotal')

MR_f_lambda_fig = figure('Color','White');
plot(vals_lambda,MR_f_lambda);
legend({'No schooling','Primary','Lower middle','Upper middle','College'}...
,'Location','northwest');
legend('boxoff')
xlabel('\lambda')
ylabel('Marriage rate per period')

MR_m_lambda_fig = figure('Color','White');
plot(vals_lambda,MR_m_lambda);
legend({'No schooling','Primary','Lower middle','Upper middle','College'}...
,'Location','northwest');
legend('boxoff')
xlabel('\lambda')
ylabel('Marriage rate per period')

marr_rates_lambda_fig = figure('Color','White');
plot(vals_lambda,[mr_f_lambda; mr_m_lambda]);
legend({'Females','Males'},'Location','northwest')
legend('boxoff')
xlabel('\lambda')
ylabel('Marriage rate per period')

match_prob_lambda_fig = figure('Color','White');
plot(vals_lambda,[pi_f_lambda; pi_m_lambda]);
legend({'Females','Males'},'Location','northwest')
legend('boxoff')
xlabel('\lambda')
ylabel('Matching probability')

exp_single_lambda_fig = figure('Color','White');
plot(vals_lambda,[exp_single_f_lambda; exp_single_m_lambda]);
legend({'Females','Males'},'Location','southwest')
legend('boxoff')
xlabel('\lambda')
ylabel('Expected number of periods spent single')

lf_lambda_fig = figure('Color','White');
plot(vals_lambda,lf_bar_lambda);
legend('boxoff')
xlabel('\lambda')
ylabel('Married female market hours as fraction of husbands')

% Changing the household upkeep parameter
param.lambda = 0.5;

vals_cbar = linspace(0.01,1,100);
id_cbar = zeros(size(vals_cbar));
delta_cbar = zeros(size(vals_cbar));
hyper_f_cbar = zeros(size(vals_cbar));
hyper_m_cbar = zeros(size(vals_cbar));
MR_f_cbar = zeros(size(Z_f,1),size(vals_cbar,2)); 
MR_m_cbar = zeros(size(Z_m,1),size(vals_cbar,2));
mr_f_cbar = zeros(size(vals_cbar));
mr_m_cbar = zeros(size(vals_cbar));
exp_single_f_cbar = zeros(size(vals_cbar));
exp_single_m_cbar = zeros(size(vals_cbar));
lf_bar_cbar = zeros(size(vals_cbar));
pi_f_cbar = zeros(size(vals_cbar));
pi_m_cbar = zeros(size(vals_cbar));

for i=1:length(vals_cbar)
    
    param.cbar = vals_cbar(i);
    [Q, ID] = marriage_eq(param,OMEGA_f,OMEGA_m,Z_f,Z_m);
    
    % Percentage of matches in which female is binding party
    id_cbar(i) = sum(sum(Z_f*Z_m'.*(ID>=0)));
    
    % Marital sorting
    MS = (Z_f*Z_m'.*(1-normcdf(Q,param.MU)))...
        /sum(sum(Z_f*Z_m'.*(1-normcdf(Q,param.MU))));
    delta_cbar(i) = sort_cont_mat(MS);
    [hyper_f_cbar(i), hyper_m_cbar(i)] = hypergamy(MS);
    
    % Marriage rates
    [pi_f,pi_m] = match_prob(param);
    MR_f_cbar(:,i) = ...
        sum(pi_f*Z_f*Z_m'.*(1-normcdf(Q,param.MU)),2)./Z_f;
    MR_m_cbar(:,i) = ...
        sum(pi_m*Z_f*Z_m'.*(1-normcdf(Q,param.MU)),1)'./Z_m;


    mr_f_cbar(i) = pi_f*Z_f'*(1-normcdf(Q,param.MU))*Z_m;
    mr_m_cbar(i) = pi_m*Z_f'*(1-normcdf(Q,param.MU))*Z_m;
    % Expected number of periods spent in single's pool
    exp_single_f_cbar(i) = 1/mr_f_cbar(i);
    exp_single_m_cbar(i) = 1/mr_m_cbar(i);
    
    % Aggregate labor supply of married women relative to men
    Lf = married_fls(param,OMEGA_f,OMEGA_m);
    lf_bar_cbar(i) = sum(sum(Lf.*MS));
    
    % Match probabilities
    [pi_f_cbar(i),pi_m_cbar(i)] = match_prob(param);
    
end

% Graphs

sorting_cbar_fig = figure('Color','White');
plot(vals_cbar,delta_cbar);
legend('boxoff')
xlabel('Household upkeep')
yticks(1)
ylabel('\delta')

hypergamy_cbar_fig = figure('Color','White');
plot(vals_cbar,[hyper_f_cbar; hyper_m_cbar])
legend({'Female', 'Male'},'Location','northeast')
legend('boxoff')
xlabel('Household upkeep')
yticks(1)
ylabel('Hypergamy measure')

id_cbar_fig = figure('Color','White');
plot(vals_cbar,id_cbar);
xlabel('Household upkeep')
ylabel('Fraction of matches female is pivotal')

MR_f_cbar_fig = figure('Color','White');
plot(vals_cbar,MR_f_cbar);
legend({'No schooling','Primary','Lower middle','Upper middle','College'}...
,'Location','southwest');
legend('boxoff')
xlabel('Household upkeep')
ylabel('Marriage rate per period')

MR_m_cbar_fig = figure('Color','White');
plot(vals_cbar,MR_m_cbar);
legend({'No schooling','Primary','Lower middle','Upper middle','College'}...
,'Location','southwest');
legend('boxoff')
xlabel('Household upkeep')
ylabel('Marriage rate per period')

marr_rates_cbar_fig = figure('Color','White');
plot(vals_cbar,[mr_f_cbar; mr_m_cbar]);
legend({'Females','Males'},'Location','southwest')
legend('boxoff')
xlabel('Household upkeep')
ylabel('Marriage rate per period')

match_prob_cbar_fig = figure('Color','White');
plot(vals_cbar,[pi_f_cbar; pi_m_cbar]);
legend({'Females','Males'},'Location','northwest')
legend('boxoff')
xlabel('Household upkeep')
ylabel('Matching probability')

exp_single_cbar_fig = figure('Color','White');
plot(vals_cbar,[exp_single_f_cbar; exp_single_m_cbar]);
legend({'Females','Males'},'Location','northwest')
legend('boxoff')
xlabel('Household upkeep')
ylabel('Expected number of periods spent single')

lf_cbar_fig = figure('Color','White');
plot(vals_cbar,lf_bar_cbar);
legend('boxoff')
xlabel('Household upkeep')
ylabel('Married female market hours as fraction of husbands')

% Changing the sex ratio among singles
param.A = 0.75;
param.lambda = 0.75;
param.cbar = 0.5;

vals_theta_S = linspace(0.4,2.5,211);
id_theta_S = zeros(size(vals_theta_S));
delta_theta_S = zeros(size(vals_theta_S));
hyper_f_theta_S = zeros(size(vals_theta_S));
hyper_m_theta_S = zeros(size(vals_theta_S));
MR_f_theta_S = zeros(size(Z_f,1),size(vals_theta_S,2)); 
MR_m_theta_S = zeros(size(Z_m,1),size(vals_theta_S,2));
mr_f_theta_S = zeros(size(vals_theta_S));
mr_m_theta_S = zeros(size(vals_theta_S));
exp_single_f_theta_S = zeros(size(vals_theta_S));
exp_single_m_theta_S = zeros(size(vals_theta_S));
lf_bar_theta_S = zeros(size(vals_theta_S));
pi_f_theta_S = zeros(size(vals_theta_S));
pi_m_theta_S = zeros(size(vals_theta_S));

for i=1:length(vals_theta_S)
    
    param.theta_S = vals_theta_S(i);
    [Q, ID] = marriage_eq(param,OMEGA_f,OMEGA_m,Z_f,Z_m);
    
    % Percentage of matches in which female is binding party
    id_theta_S(i) = sum(sum(Z_f*Z_m'.*(ID>=0)));
    
    % Marital sorting
    MS = (Z_f*Z_m'.*(1-normcdf(Q,param.MU)))...
        /sum(sum(Z_f*Z_m'.*(1-normcdf(Q,param.MU))));
    delta_theta_S(i) = sort_cont_mat(MS);
    [hyper_f_theta_S(i), hyper_m_theta_S(i)] = hypergamy(MS);
    
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

% Graphs

sorting_theta_S_fig = figure('Color','White');
plot(vals_theta_S,delta_theta_S);
legend('boxoff')
xlabel('\theta_S')
yticks(1)
ylabel('\delta')

hypergamy_theta_S_fig = figure('Color','White');
plot(vals_theta_S,[hyper_f_theta_S; hyper_m_theta_S])
legend({'Female', 'Male'},'Location','northeast')
legend('boxoff')
xlabel('\theta_S')
yticks(1)
ylabel('Hypergamy measure')

id_theta_S_fig = figure('Color','White');
plot(vals_theta_S,id_theta_S);
xlabel('\theta_S')
ylabel('Fraction of matches female is pivotal')

MR_f_theta_S_fig = figure('Color','White');
plot(vals_theta_S,MR_f_theta_S);
legend({'No schooling','Primary','Lower middle','Upper middle','College'}...
,'Location','southwest');
legend('boxoff')
xlabel('\theta_S')
ylabel('Marriage rate per period')

MR_m_theta_S_fig = figure('Color','White');
plot(vals_theta_S,MR_m_theta_S);
legend({'No schooling','Primary','Lower middle','Upper middle','College'}...
,'Location','southwest');
legend('boxoff')
xlabel('\theta_S')
ylabel('Marriage rate per period')

marr_rates_theta_S_fig = figure('Color','White');
plot(vals_theta_S,[mr_f_theta_S; mr_m_theta_S]);
legend({'Females','Males'},'Location','southwest')
legend('boxoff')
xlabel('\theta_S')
ylabel('Marriage rate per period')

match_prob_theta_S_fig = figure('Color','White');
plot(vals_theta_S,[pi_f_theta_S; pi_m_theta_S]);
legend({'Females','Males'},'Location','northwest')
legend('boxoff')
xlabel('\theta_S')
ylabel('Matching probability')

exp_single_theta_S_fig = figure('Color','White');
plot(vals_theta_S,[exp_single_f_theta_S; exp_single_m_theta_S]);
legend({'Females','Males'},'Location','northwest')
legend('boxoff')
xlabel('\theta_S')
ylabel('Expected number of periods spent single')

lf_theta_S_fig = figure('Color','White');
plot(vals_theta_S,lf_bar_theta_S);
legend('boxoff')
xlabel('\theta_S')
ylabel('Married female market hours as fraction of husbands')

% Changing the gender wage gap
param.theta_S = 1;

vals_gender_gap = linspace(0.5,1,51);
id_gender_gap = zeros(size(vals_gender_gap));
delta_gender_gap = zeros(size(vals_gender_gap));
hyper_f_gender_gap = zeros(size(vals_gender_gap));
hyper_m_gender_gap = zeros(size(vals_gender_gap));
MR_f_gender_gap = zeros(size(Z_f,1),size(vals_gender_gap,2)); 
MR_m_gender_gap = zeros(size(Z_m,1),size(vals_gender_gap,2));
mr_f_gender_gap = zeros(size(vals_gender_gap));
mr_m_gender_gap = zeros(size(vals_gender_gap));
exp_single_f_gender_gap = zeros(size(vals_gender_gap));
exp_single_m_gender_gap = zeros(size(vals_gender_gap));
lf_bar_gender_gap = zeros(size(vals_gender_gap));
pi_f_gender_gap = zeros(size(vals_gender_gap));
pi_m_gender_gap = zeros(size(vals_gender_gap));

for i=1:length(vals_gender_gap)
    
    param.gender_wage_gap = vals_gender_gap(i);
    OMEGA_f = param.gender_wage_gap*OMEGA_m;
    [Q, ID] = marriage_eq(param,OMEGA_f,OMEGA_m,Z_f,Z_m);
    
    % Percentage of matches in which female is binding party
    id_gender_gap(i) = sum(sum(Z_f*Z_m'.*(ID>=0)));
    
    % Marital sorting
    MS = (Z_f*Z_m'.*(1-normcdf(Q,param.MU)))...
        /sum(sum(Z_f*Z_m'.*(1-normcdf(Q,param.MU))));
    delta_gender_gap(i) = sort_cont_mat(MS);
    [hyper_f_gender_gap(i), hyper_m_gender_gap(i)] = hypergamy(MS);
    
    % Marriage rates
    [pi_f,pi_m] = match_prob(param);
    MR_f_gender_gap(:,i) = ...
        sum(pi_f*Z_f*Z_m'.*(1-normcdf(Q,param.MU)),2)./Z_f;
    MR_m_gender_gap(:,i) = ...
        sum(pi_m*Z_f*Z_m'.*(1-normcdf(Q,param.MU)),1)'./Z_m;


    mr_f_gender_gap(i) = pi_f*Z_f'*(1-normcdf(Q,param.MU))*Z_m;
    mr_m_gender_gap(i) = pi_m*Z_f'*(1-normcdf(Q,param.MU))*Z_m;
    % Expected number of periods spent in single's pool
    exp_single_f_gender_gap(i) = 1/mr_f_gender_gap(i);
    exp_single_m_gender_gap(i) = 1/mr_m_gender_gap(i);
    
    % Aggregate labor supply of married women relative to men
    Lf = married_fls(param,OMEGA_f,OMEGA_m);
    lf_bar_gender_gap(i) = sum(sum(Lf.*MS));
    
    % Match probabilities
    [pi_f_gender_gap(i),pi_m_gender_gap(i)] = match_prob(param);
    
end

% Graphs

sorting_gender_gap_fig = figure('Color','White');
plot(vals_gender_gap,delta_gender_gap);
legend('boxoff')
xlabel('\phi')
yticks(1)
ylabel('\delta')

hypergamy_gender_gap_fig = figure('Color','White');
plot(vals_gender_gap,[hyper_f_gender_gap; hyper_m_gender_gap])
legend({'Female', 'Male'},'Location','northeast')
legend('boxoff')
xlabel('\phi')
yticks(1)
ylabel('Hypergamy measure')

id_gender_gap_fig = figure('Color','White');
plot(vals_gender_gap,id_gender_gap);
xlabel('\phi')
ylabel('Fraction of matches female is pivotal')

MR_f_gender_gap_fig = figure('Color','White');
plot(vals_gender_gap,MR_f_gender_gap);
legend({'No schooling','Primary','Lower middle','Upper middle','College'}...
,'Location','southwest');
legend('boxoff')
xlabel('\phi')
ylabel('Marriage rate per period')

MR_m_gender_gap_fig = figure('Color','White');
plot(vals_gender_gap,MR_m_gender_gap);
legend({'No schooling','Primary','Lower middle','Upper middle','College'}...
,'Location','southwest');
legend('boxoff')
xlabel('\phi')
ylabel('Marriage rate per period')

marr_rates_gender_gap_fig = figure('Color','White');
plot(vals_gender_gap,[mr_f_gender_gap; mr_m_gender_gap]);
legend({'Females','Males'},'Location','southwest')
legend('boxoff')
xlabel('\phi')
ylabel('Marriage rate per period')

match_prob_gender_gap_fig = figure('Color','White');
plot(vals_gender_gap,[pi_f_gender_gap; pi_m_gender_gap]);
legend({'Females','Males'},'Location','northwest')
legend('boxoff')
xlabel('\phi')
ylabel('Matching probability')

exp_single_gender_gap_fig = figure('Color','White');
plot(vals_gender_gap,[exp_single_f_gender_gap; exp_single_m_gender_gap]);
legend({'Females','Males'},'Location','northwest')
legend('boxoff')
xlabel('\phi')
ylabel('Expected number of periods spent single')

lf_gender_gap_fig = figure('Color','White');
plot(vals_gender_gap,lf_bar_gender_gap);
legend('boxoff')
xlabel('\phi')
ylabel('Married female market hours as fraction of husbands')

% Adding taste for homogamy
param.gender_wage_gap = 1;
OMEGA_f = param.gender_wage_gap*OMEGA_m;

vals_taste_homogamy = linspace(0,1,101);
id_taste_homogamy = zeros(size(vals_taste_homogamy));
delta_taste_homogamy = zeros(size(vals_taste_homogamy));
hyper_f_taste_homogamy = zeros(size(vals_taste_homogamy));
hyper_m_taste_homogamy = zeros(size(vals_taste_homogamy));
MR_f_taste_homogamy = zeros(size(Z_f,1),size(vals_taste_homogamy,2)); 
MR_m_taste_homogamy = zeros(size(Z_m,1),size(vals_taste_homogamy,2));
mr_f_taste_homogamy = zeros(size(vals_taste_homogamy));
mr_m_taste_homogamy = zeros(size(vals_taste_homogamy));
exp_single_f_taste_homogamy = zeros(size(vals_taste_homogamy));
exp_single_m_taste_homogamy = zeros(size(vals_taste_homogamy));
lf_bar_taste_homogamy = zeros(size(vals_taste_homogamy));
pi_f_taste_homogamy = zeros(size(vals_taste_homogamy));
pi_m_taste_homogamy = zeros(size(vals_taste_homogamy));

for i=1:length(vals_taste_homogamy)
    
    taste_homogamy = vals_taste_homogamy(i) ;
    param.MU = taste_homogamy*eye(types_f);
    
    OMEGA_f = param.gender_wage_gap*OMEGA_m;
    [Q, ID] = marriage_eq(param,OMEGA_f,OMEGA_m,Z_f,Z_m);
    
    % Percentage of matches in which female is binding party
    id_taste_homogamy(i) = sum(sum(Z_f*Z_m'.*(ID>=0)));
    
    % Marital sorting
    MS = (Z_f*Z_m'.*(1-normcdf(Q,param.MU)))...
        /sum(sum(Z_f*Z_m'.*(1-normcdf(Q,param.MU))));
    delta_taste_homogamy(i) = sort_cont_mat(MS);
    [hyper_f_taste_homogamy(i), hyper_m_taste_homogamy(i)] = hypergamy(MS);
    
    % Marriage rates
    [pi_f,pi_m] = match_prob(param);
    MR_f_taste_homogamy(:,i) = ...
        sum(pi_f*Z_f*Z_m'.*(1-normcdf(Q,param.MU)),2)./Z_f;
    MR_m_taste_homogamy(:,i) = ...
        sum(pi_m*Z_f*Z_m'.*(1-normcdf(Q,param.MU)),1)'./Z_m;


    mr_f_taste_homogamy(i) = pi_f*Z_f'*(1-normcdf(Q,param.MU))*Z_m;
    mr_m_taste_homogamy(i) = pi_m*Z_f'*(1-normcdf(Q,param.MU))*Z_m;
    % Expected number of periods spent in single's pool
    exp_single_f_taste_homogamy(i) = 1/mr_f_taste_homogamy(i);
    exp_single_m_taste_homogamy(i) = 1/mr_m_taste_homogamy(i);
    
    % Aggregate labor supply of married women relative to men
    Lf = married_fls(param,OMEGA_f,OMEGA_m);
    lf_bar_taste_homogamy(i) = sum(sum(Lf.*MS));
    
    % Match probabilities
    [pi_f_taste_homogamy(i),pi_m_taste_homogamy(i)] = match_prob(param);
    
end

% Graphs

sorting_taste_homogamy_fig = figure('Color','White');
plot(vals_taste_homogamy,delta_taste_homogamy);
legend('boxoff')
xlabel('\kappa')
ylabel('\delta')

hypergamy_taste_homogamy_fig = figure('Color','White');
plot(vals_taste_homogamy,[hyper_f_taste_homogamy; hyper_m_taste_homogamy])
legend({'Female', 'Male'},'Location','northeast')
legend('boxoff')
xlabel('\kappa')
ylabel('Hypergamy measure')

id_taste_homogamy_fig = figure('Color','White');
plot(vals_taste_homogamy,id_taste_homogamy);
xlabel('\kappa')
ylabel('Fraction of matches female is pivotal')

MR_f_taste_homogamy_fig = figure('Color','White');
plot(vals_taste_homogamy,MR_f_taste_homogamy);
legend({'No schooling','Primary','Lower middle','Upper middle','College'}...
,'Location','southwest');
legend('boxoff')
xlabel('\kappa')
ylabel('Marriage rate per period')

MR_m_taste_homogamy_fig = figure('Color','White');
plot(vals_taste_homogamy,MR_m_taste_homogamy);
legend({'No schooling','Primary','Lower middle','Upper middle','College'}...
,'Location','southwest');
legend('boxoff')
xlabel('\kappa')
ylabel('Marriage rate per period')

marr_rates_taste_homogamy_fig = figure('Color','White');
plot(vals_taste_homogamy,[mr_f_taste_homogamy; mr_m_taste_homogamy]);
legend({'Females','Males'},'Location','southwest')
legend('boxoff')
xlabel('\kappa')
ylabel('Marriage rate per period')

match_prob_taste_homogamy_fig = figure('Color','White');
plot(vals_taste_homogamy,[pi_f_taste_homogamy; pi_m_taste_homogamy]);
legend({'Females','Males'},'Location','northwest')
legend('boxoff')
xlabel('\kappa')
ylabel('Matching probability')

exp_single_taste_homogamy_fig = figure('Color','White');
plot(vals_taste_homogamy,[exp_single_f_taste_homogamy; ...
    exp_single_m_taste_homogamy]);
legend({'Females','Males'},'Location','northwest')
legend('boxoff')
xlabel('\kappa')
ylabel('Expected number of periods spent single')

lf_taste_homogamy_fig = figure('Color','White');
plot(vals_taste_homogamy,lf_bar_taste_homogamy);
legend('boxoff')
xlabel('\kappa')
ylabel('Married female market hours as fraction of husbands')

%{
Scenario 2: income inequality
%}



% % Baseline paramaters
% 
% periods_year = 1;
% param.beta = (0.96)^(1/periods_year); % Discount rate
% life_expectancy = 45*periods_year; % Life expectancy at age 20
% param.delta = 1/life_expectancy; % Death rate
% param.A = 0.7; % Efficiency of matching function
% param.alpha = 0.5; % Elasticity of matching function
% param.MU = 0.5*eye(5); % Mean of the distribution of match quality by 
%                        % partners' education 
% param.cbar = 0.75; % Household upkeep cost
% param.lambda = 0.75; % Weight of market goods consumption in household 
%                   % utility
% param.theta_S = 1.34; % Sex ratio among singles
% param.skill_premium = linspace(1,2,5)'; % Skill premium
% param.gender_wage_gap = 0.75; % Gender wage gap
% 
% % From the data
% 
% % Distribution of education among singles
% Z_f = [0.0932 ; 0.087 ; 0.3478 ; 0.2981 ; 0.1739];
% Z_m = [0.1019 ; 0.088 ; 0.4167 ; 0.2730 ; 0.1204];
% 
% % Z_f = ones(5,1)/5;
% % Z_m = Z_f;
% 
% % Skill premium vectors
% OMEGA_m = param.skill_premium;
% OMEGA_f = param.gender_wage_gap*OMEGA_m;
% 
% % Changing the sex ratio among singles
% values_theta_S = linspace(0.5,2.5,201);
% id_theta_S = zeros(size(values_theta_S));
% delta_theta_S = zeros(size(values_theta_S));
% hyper_f_theta_S = zeros(size(values_theta_S));
% hyper_m_theta_S = zeros(size(values_theta_S));
% MR_f_theta_S = zeros(size(Z_f,1),size(values_theta_S,2)); 
% MR_m_theta_S = zeros(size(Z_m,1),size(values_theta_S,2));
% mr_f_theta_S = zeros(size(values_theta_S));
% mr_m_theta_S = zeros(size(values_theta_S));
% exp_single_f_theta_S = zeros(size(values_theta_S));
% exp_single_m_theta_S = zeros(size(values_theta_S));
% lf_bar_theta_S = zeros(size(values_theta_S));
% pi_f_theta_S = zeros(size(values_theta_S));
% pi_m_theta_S = zeros(size(values_theta_S));
% 
% for i=1:length(values_theta_S)
%     
%     param.theta_S = values_theta_S(i);
%     [Q, ID] = marriage_eq(param,OMEGA_f,OMEGA_m,Z_f,Z_m);
%     
%     % Percentage of matches in which female is binding party
%     id_theta_S(i) = sum(sum(Z_f*Z_m'.*ID));
%     
%     % Marital sorting
%     MS = (Z_f*Z_m'.*(1-normcdf(Q,param.MU)))...
%         /sum(sum(Z_f*Z_m'.*(1-normcdf(Q,param.MU))));
%     delta_theta_S(i) = sort_cont_mat(MS);
%     [hyper_f_theta_S(i), hyper_m_theta_S(i)] = hypergamy(MS);
%     
%     % Marriage rates
%     [pi_f,pi_m] = match_prob(param);
%     MR_f_theta_S(:,i) = ...
%         sum(pi_f*Z_f*Z_m'.*(1-normcdf(Q,param.MU)),2)./Z_f;
%     MR_m_theta_S(:,i) = ...
%         sum(pi_m*Z_f*Z_m'.*(1-normcdf(Q,param.MU)),1)'./Z_m;
% 
% 
%     mr_f_theta_S(i) = pi_f*Z_f'*(1-normcdf(Q,param.MU))*Z_m;
%     mr_m_theta_S(i) = pi_m*Z_f'*(1-normcdf(Q,param.MU))*Z_m;
%     % Expected number of periods spent in single's pool
%     exp_single_f_theta_S(i) = 1/mr_f_theta_S(i);
%     exp_single_m_theta_S(i) = 1/mr_m_theta_S(i);
%     
%     % Aggregate labor supply of married women relative to men
%     Lf = married_fls(param,OMEGA_f,OMEGA_m);
%     lf_bar_theta_S(i) = sum(sum(Lf.*MS));
%     
%     % Match probabilities
%     [pi_f_theta_S(i),pi_m_theta_S(i)] = match_prob(param);
%     
% end
% 
% sorting_theta_S_fig = figure('Color','White');
% plot(values_theta_S,delta_theta_S);
% legend('boxoff')
% xlabel('\theta_S')
% 
% hypergamy_theta_S_fig = figure('Color','White');
% plot(values_theta_S,[hyper_f_theta_S; hyper_m_theta_S])
% legend({'Female', 'Male'},'Location','northeast')
% legend('boxoff')
% xlabel('\theta_S')
% 
% id_theta_S_fig = figure('Color','White');
% plot(values_theta_S,id_theta_S);
% xlabel('\theta_S')
% 
% MR_f_theta_S_fig = figure('Color','White');
% plot(values_theta_S,MR_f_theta_S);
% legend({'No schooling','Primary','Lower middle','Upper middle','College'}...
% ,'Location','northwest');
% legend('boxoff')
% xlabel('\theta_S')
% 
% MR_m_theta_S_fig = figure('Color','White');
% plot(values_theta_S,MR_m_theta_S);
% legend({'No schooling','Primary','Lower middle','Upper middle','College'}...
% ,'Location','northeast');
% legend('boxoff')
% xlabel('\theta_S')
% 
% marr_rates_theta_S_fig = figure('Color','White');
% plot(values_theta_S,[mr_f_theta_S; mr_m_theta_S]);
% legend({'Females','Males'},'Location','northwest')
% legend('boxoff')
% xlabel('\theta_S')
% 
% match_prob_theta_S_fig = figure('Color','White');
% plot(values_theta_S,[pi_f_theta_S; pi_m_theta_S]);
% legend({'Females','Males'},'Location','northwest')
% legend('boxoff')
% xlabel('\theta_S')
% 
% lf_theta_S_fig = figure('Color','White');
% plot(values_theta_S,lf_bar_theta_S);
% legend('boxoff')
% xlabel('\theta_S')


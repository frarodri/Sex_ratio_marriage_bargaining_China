%**************************************************************************
% Comparative statics - two parameters                                    *
% Francisco Javier Rodríguez Román                                        * 
% Department of Economics                                                 *
% Universidad Carlos III de Madrid                                        *
%**************************************************************************

% home: G:\My Drive\China project\Model
% uni: G:\Mi unidad\China project\Model

clear ; close all; clc
cd('G:\Mi unidad\China project\Model')

% Choose settings for the comparative statics exercise

types = {'Skill levels','Percentiles in hourly wage distribution'};

[agent_types,~] = listdlg('PromptString','Agent types correspond to:',...
    'ListString',types,'SelectionMode','single'); 

list = {'A parameter','Expectations','The wage structure'};

[type_first,~] = listdlg('PromptString',...
                         'Perform comparative statics primarily on:',...
                         'ListString',list,'SelectionMode','single');

promt = 'Perform comparative statics primarily specifically on: ';
first = input(promt,'s');

promt2 = 'Minimum value: ';
min_first = input(promt2);

promt3 = 'Step: ';
step_first = input(promt3);

promt4 = 'Maximum value: ';
max_first = input(promt4);

[type_second,~] = listdlg('PromptString',...
                          'Perform comparative statics second on:',...
                          'ListString',list,'SelectionMode','single');

promt5 = 'Perform comparative statics second specifically on: ';
second = input(promt5,'s');

promt6 = 'Vector of values for second: ';
vals2 = input(promt6,'s');

promt7 = 'Title for x-axis in graphs: ';
title_x_axis = input(promt7,'s');

promt8 = 'Variable name in legend: ';
var_legend = input(promt8,'s');

save_graphs = questdlg('Save graphs?','Save graphs?');

if strcmp(save_graphs,'Yes')==1
    
    promt9 = 'Graph suffix: ';
    graph_suffix = input(promt9,'s');
    
end

% Baseline parameters
%%%%%%%%%%%%%%%%%%%%%

if agent_types == 1
    
    types_number = 5;
    
elseif agent_types ==2
    
    types_number = 100;
    
end

param.periods_year = 1;
param.beta = (0.96)^(1/param.periods_year); % Discount rate

life_expectancy = 45*param.periods_year; % Life expectancy at age 20
param.delta = 1/life_expectancy; % Death rate

param.A = 0.8; % Efficiency of matching function
param.alpha = 0.5; % Elasticity of matching function

param.kappa = 0;
param.MU = param.kappa*eye(types_number); % Mean of the distribution of 
                                          % match quality by partners' 
                                          % education 

param.cbar = 0.5; % Household upkeep cost/minimum survival consumption
param.lambda = 0.8; % Weight of market goods consumption in household 
                    % utility
param.gamma = 1.5; % Economies of scale in minimum survival consumption
                   % (1 is perfect pubic good, 2 is no economies of scale) 
                   
% Baseline expectations
%%%%%%%%%%%%%%%%%%%%%%%

% Sex ratio among singles
expect.theta_S = 1.5; 

% Distribution of types among singles

expect.P_f = ones(types_number,1)/types_number; 
expect.P_m = ones(types_number,1)/types_number;

% Baseline wage structure
%%%%%%%%%%%%%%%%%%%%%%%%%

wages.gender_ratio = 0.75;
wages.skill_premium = 0.5;
wages.general_level = 1;
wages.skills = (1:1:types_number)';
wages.k = 1.1;
wages.sigma = 0.5;
wages.theta = 1;

if agent_types == 1
    
    wages.OMEGA_m = wages.general_level*...
        (ones(size(expect.P_m))+wages.skill_premium*(wages.skills-1));
    wages.OMEGA_f = wages.gender_ratio*wages.OMEGA_m;
        
elseif agent_types == 2
    
    wages.OMEGA_m = ...
        gpinv(linspace(0,0.99,100),wages.k,wages.sigma,wages.theta)';
    wages.OMEGA_f = wages.gender_ratio*wages.OMEGA_m;
    
end

% Pre-allocate vectors
%%%%%%%%%%%%%%%%%%%%%%

values_first = (min_first:step_first:max_first);
values_second = str2num(vals2);

VS_f_vec = zeros(size(expect.P_f,1),size(values_first,2),...
    size(values_second,2));

VS_m_vec = zeros(size(expect.P_m,1),size(values_first,2),...
    size(values_second,2));

id = zeros(size(values_first,2), size(values_second,2));

delta_vec = zeros(size(values_first,2), size(values_second,2));
hypgam_f_vec = zeros(size(values_first,2), size(values_second,2));
hypgam_m_vec = zeros(size(values_first,2), size(values_second,2));

delta_2_vec = zeros(size(values_first,2), size(values_second,2));
hypgam_f_2_vec = zeros(size(values_first,2), size(values_second,2));
hypgam_m_2_vec = zeros(size(values_first,2), size(values_second,2));

MR_f = zeros(size(expect.P_f,1),size(values_first,2),...
    size(values_second,2)); 
MR_m = zeros(size(expect.P_m,1),size(values_first,2),...
    size(values_second,2));

mr_f = zeros(size(values_first,2), size(values_second,2));
mr_m = zeros(size(values_first,2), size(values_second,2));

exp_single_f = zeros(size(values_first,2), size(values_second,2));
exp_single_m = zeros(size(values_first,2), size(values_second,2));

lf_bar = zeros(size(values_first,2), size(values_second,2));

pi_f_vec = zeros(size(values_first,2), size(values_second,2));
pi_m_vec = zeros(size(values_first,2), size(values_second,2));

ineq_inc = zeros(size(values_first,2), size(values_second,2));
ineq_cons = zeros(size(values_first,2), size(values_second,2));

for i = 1:size(values_first,2)

    if type_first == 1
              
        param.(first) = values_first(i);
        param.MU = param.kappa*eye(types_number);
        
    elseif type_first == 2
        
        expect.(first) = values_first(i);
        
    elseif type_first == 3
        
        wages.(first) = values_first(i);
        
        if agent_types == 1
            
            wages.OMEGA_m = wages.general_level*...
                (ones(size(expect.P_m))+...
                wages.skill_premium*(wages.skills-1));
            wages.OMEGA_f = wages.gender_ratio*wages.OMEGA_m;

        elseif agent_types == 2

            wages.OMEGA_m = gpinv(linspace(0,0.99,100),...
                wages.k,wages.sigma,wages.theta)';
            wages.OMEGA_f = wages.gender_ratio*wages.OMEGA_m;

        end
        
    end
    
    for j = 1:size(values_second,2)
        
       if type_second == 1

            param.(second) = values_second(j);
            param.MU = param.kappa*eye(types_number);

       elseif type_second == 2

            expect.(second) = values_second(j);

       elseif type_second == 3

            wages.(second) = values_second(j);

            if agent_types == 1

                wages.OMEGA_m = ...
                    wages.general_level*(ones(size(expect.P_m))+...
                    wages.skill_premium*(wages.skills-1));
                wages.OMEGA_f = wages.gender_ratio*wages.OMEGA_m;

            elseif agent_types == 2

                wages.OMEGA_m = gpinv(linspace(0,0.99,100),...
                                      wages.k,wages.sigma,wages.theta)';
                wages.OMEGA_f = wages.gender_ratio*wages.OMEGA_m;

            end

        end
       
        [Q, VS_f, VS_m, ID, MS] = marriage_pree(param,expect,wages);
        
        % Value of being single
        VS_f_vec(:,i,j) = VS_f;
        VS_m_vec(:,i,j) = VS_m;
        
        % Marital sorting
        [delta, h_f, h_m] = sort_cont_mat(MS);
        delta_vec(i,j) = delta;
        hypgam_f_vec(i,j) = h_f;
        hypgam_m_vec(i,j) = h_m;
        
        [delta2, h_f2, h_m2] = sort_cont_mat_2(MS,expect.P_f,expect.P_m);
        delta_2_vec(i,j) = delta2;
        hypgam_f_2_vec(i,j) = h_f2;
        hypgam_m_2_vec(i,j) = h_m2;
        
        % Match probabilities
        [pi_f,pi_m] = match_prob(param,expect.theta_S);
        pi_f_vec(i,j) = pi_f;
        pi_m_vec(i,j) = pi_m;

        % Marriage rates
        MR_f(:,i,j) = pi_f*(1-normcdf(Q,param.MU))*expect.P_m;
        MR_m(:,i,j) = pi_m*(1-normcdf(Q,param.MU))'*expect.P_f;

        mr_f(i,j) = expect.P_f'*MR_f(:,i,j);
        mr_m(i,j) = expect.P_m'*MR_m(:,i,j);

        % Expected number of periods spent in single's pool
        exp_single_f(i,j) = 1/mr_f(i,j);
        exp_single_m(i,j) = 1/mr_m(i,j);

        % Labor supply of married women relative to men
        Lf = married_fls(param,wages);
        lf_bar(i,j) = sum(sum(Lf.*MS));

        % Household inequality
        [ineq_inc(i,j), ineq_cons(i,j)] = ...
        ineq(param,wages,Q,MS,expect.theta_S,expect.P_f,expect.P_m);
       
    end
    
end
   
legend_text = string(var_legend)+"="+values_second;

match_prob_f_fig = figure('Color','White');
plot(values_first,pi_f_vec);
legend(legend_text,'Interpreter','latex','Location','best')
legend('boxoff')
xlabel(title_x_axis,'interpreter','latex')
ylabel('Matching probability')

match_prob_m_fig = figure('Color','White');
plot(values_first,pi_m_vec);
legend(legend_text,'Interpreter','latex','Location','best')
legend('boxoff')
xlabel(title_x_axis,'interpreter','latex')
ylabel('Matching probability')

marr_rates_f_fig = figure('Color','White');
plot(values_first,mr_f);
legend(legend_text,'Interpreter','latex','Location','best')
legend('boxoff')
xlabel(title_x_axis,'interpreter','latex')
ylabel('Marriage rate per period')

marr_rates_m_fig = figure('Color','White');
plot(values_first,mr_m);
legend(legend_text,'Interpreter','latex','Location','best')
legend('boxoff')
xlabel(title_x_axis,'interpreter','latex')
ylabel('Marriage rate per period')

exp_single_f_fig = figure('Color','White');
plot(values_first,exp_single_f);
legend(legend_text,'Interpreter','latex','Location','best')
legend('boxoff')
xlabel(title_x_axis,'interpreter','latex')
ylabel('Expected number of periods spent single')

exp_single_m_fig = figure('Color','White');
plot(values_first,exp_single_m);
legend(legend_text,'Interpreter','latex','Location','best')
legend('boxoff')
xlabel(title_x_axis,'interpreter','latex')
ylabel('Expected number of periods spent single')

sorting_fig = figure('Color','White');
plot(values_first,delta_vec);
legend(legend_text,'Interpreter','latex','Location','best')
legend('boxoff')
xlabel(title_x_axis,'interpreter','latex')
ylabel('\delta')

hypergamy_f_fig = figure('Color','White');
plot(values_first,hypgam_f_vec)
legend(legend_text,'Interpreter','latex','Location','best')
legend('boxoff')
xlabel(title_x_axis,'interpreter','latex')
ylabel('Hypergamy measure')

hypergamy_m_fig = figure('Color','White');
plot(values_first,hypgam_m_vec)
legend(legend_text,'Interpreter','latex','Location','best')
legend('boxoff')
xlabel(title_x_axis,'interpreter','latex')
ylabel('Hypergamy measure')

lf_fig = figure('Color','White');
plot(values_first,lf_bar);
legend(legend_text,'Interpreter','latex','Location','best')
legend('boxoff')
xlabel(title_x_axis,'interpreter','latex')
ylabel('Married female market hours as fraction of husbands')

ineq_fig = figure('Color','White');
plot(values_first,ineq_inc);
legend(legend_text,'Interpreter','latex','Location','best')
legend('boxoff')
xlabel(title_x_axis,'interpreter','latex')
ylabel('Gini index')


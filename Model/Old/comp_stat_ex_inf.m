%**************************************************************************
% Comparative statics - exogenous inflow of singles                       *
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
[index,~] = listdlg('PromptString','Perform comparative statics on:',...
                    'ListString',list,'SelectionMode','single');

promt = 'Perform comparative statics specifically on: ';
changing = input(promt,'s');

promt2 = 'Minimum value: ';
min = input(promt2);

promt3 = 'Step: ';
step = input(promt3);

promt4 = 'Maximum value: ';
max = input(promt4);

promt5 = 'Title for x-axis in graphs: ';
title_x_axis = input(promt5,'s');

save_graphs = questdlg('Save graphs?','Save graphs?');

if strcmp(save_graphs,'Yes')==1
    
    promt6 = 'Graph suffix: ';
    graph_suffix = input(promt6,'s');
    
end

% Baseline parameters
%%%%%%%%%%%%%%%%%%%%%

if agent_types == 1
    
    types_number = 5;
    
elseif agent_types == 2
    
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
                 
param.theta = 1.2; % Sex ratio among new entrants/general population
param.P_f = ones(types_number,1)/types_number; % Distribution over types
                                               % among new female entrants                                               
param.P_m = ones(types_number,1)/types_number; % Distribution over types
                                               % among new male entrants
                   
% Baseline expectations
%%%%%%%%%%%%%%%%%%%%%%%

% Sex ratio among singles
expect.theta_S = param.theta; 

% Distribution of types among singles
expect.P_f = param.P_f; 
expect.P_m = param.P_m;

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
        (ones(size(param.P_m))+wages.skill_premium*(wages.skills-1));
    wages.OMEGA_f = wages.gender_ratio*wages.OMEGA_m;
        
elseif agent_types == 2
    
    wages.OMEGA_m = ...
        gpinv(linspace(0,0.99,100),wages.k,wages.sigma,wages.theta)';
    wages.OMEGA_f = wages.gender_ratio*wages.OMEGA_m;
    
end

% Pre-allocate vectors
%%%%%%%%%%%%%%%%%%%%%%

values = (min:step:max);

theta_S_vec = zeros(size(values));
P_f_S_vec = zeros(size(param.P_f,1),size(values,2));
P_m_S_vec = zeros(size(param.P_m,1),size(values,2));

VS_f_vec = zeros(length(expect.P_f),length(values));
VS_m_vec = zeros(length(expect.P_m),length(values));

id = zeros(size(values));

delta_vec = zeros(size(values));
hypgam_f_vec = zeros(size(values));
hypgam_m_vec = zeros(size(values));

delta_2_vec = zeros(size(values));
hypgam_f_2_vec = zeros(size(values));
hypgam_m_2_vec = zeros(size(values));

MR_f = zeros(length(expect.P_f),length(values));
MR_m = zeros(length(expect.P_m),length(values));
mr_f = zeros(size(values));
mr_m = zeros(size(values));
exp_single_f = zeros(size(values));
exp_single_m = zeros(size(values));

lf_bar = zeros(size(values));

pi_f_vec = zeros(size(values));
pi_m_vec = zeros(size(values));

ineq_inc = zeros(size(values));
ineq_cons = zeros(size(values));

for i = 1:length(values)
     
  if index==1
      
      param.(changing) = values(i);
      param.MU = param.kappa*eye(types_number);
            
      expect.theta_S = param.theta; 
      expect.P_f = param.P_f; 
      expect.P_m = param.P_m;
      
  elseif index==2
      
      expect.(changing) = values(i);

  elseif index == 3
        
        wages.(first) = values_first(i);
        
        if agent_types == 1
    
            wages.OMEGA_m = wages.general_level*(ones(size(param.P_m))+...
                            wages.skill_premium*(wages.skills-1));
            wages.OMEGA_f = wages.gender_gap*wages.OMEGA_m;

        elseif agent_types == 2

            wages.OMEGA_m = gpinv(linspace(0,0.99,100),...
                                  wages.k,wages.sigma,wages.theta)';
            wages.OMEGA_f = wages.gender_gap*wages.OMEGA_m;

        end
        
   end
  
   % Compute PREE 
   [Q, ~, ~, ID, ~] = marriage_pree(param,expect,wages);
   
   % Compute the actual stationary sex ratio and distribution of types over 
   % singles
   [theta_S, P_f, P_m] = singles_dist(Q,param);
   theta_S_vec(i) = theta_S;
   P_f_S_vec(:,i) = P_f;
   P_m_S_vec(:,i) = P_m;
   
   % Compute the actual value of being single
   
   
   % Compute the actual marital sorting matrix
   MS = (P_f*P_m').*(1-normcdf(Q,param.MU))/...
            sum(sum((P_f*P_m').*(1-normcdf(Q,param.MU))));
        
   % Marital sorting
   [delta, h_f, h_m] = sort_cont_mat(MS);
   delta_vec(i) = delta;
   hypgam_f_vec(i) = h_f;
   hypgam_m_vec(i) = h_m;

   [delta2, h_f2, h_m2] = sort_cont_mat_2(MS,expect.P_f,expect.P_m);
   delta_2_vec(i) = delta2;
   hypgam_f_2_vec(i) = h_f2;
   hypgam_m_2_vec(i) = h_m2;
  
   % Match probabilities
   [pi_f,pi_m] = match_prob(param,theta_S);
   pi_f_vec(i) = pi_f;
   pi_m_vec(i) = pi_m;
  
   % Marriage rates
   MR_f(:,i) = pi_f*(1-normcdf(Q,param.MU))*P_m;
   MR_m(:,i) = pi_m*(1-normcdf(Q,param.MU))'*P_f;
   
   mr_f(i) = P_f'*MR_f(:,i);
   mr_m(i) = P_m'*MR_m(:,i);
  
   % Expected number of periods spent in single's pool
   exp_single_f(i) = 1/mr_f(i);
   exp_single_m(i) = 1/mr_m(i);
  
   % Labor supply of married women relative to men
   Lf = married_fls(param,wages);
   lf_bar(i) = sum(sum(Lf.*MS));
  
   % Household inequality
   [ineq_inc(i), ineq_cons(i)] = ...
     ineq(param,wages,Q,MS,theta_S,P_f,P_m);
  
end

% VS_f_fig = figure('Color','White');
% plot(values,VS_f_vec);
% legend({'No schooling','Primary','Lower middle','Upper middle','College'}...
% ,'Location','best');
% legend('boxoff')
% xlabel(title_x_axis,'interpreter','latex')
% ylabel('Value of being single')
% 
% VS_m_fig = figure('Color','White');
% plot(values,VS_m_vec);
% legend({'No schooling','Primary','Lower middle','Upper middle','College'}...
% ,'Location','best');
% legend('boxoff')
% xlabel(title_x_axis,'interpreter','latex')
% ylabel('Value of being single')

theta_S_fig = figure('Color','White');
plot(values,theta_S_vec);
xlabel(title_x_axis,'interpreter','latex')
ylabel('\theta_S')

if agent_types == 1
    
    P_f_S_fig = figure('Color','White');
    plot(values,P_f_S_vec);
    legend({'No schooling','Primary','Lower middle','Upper middle',...
            'College'},'Location','best');
    legend('boxoff')
    xlabel(title_x_axis,'interpreter','latex')
    ylabel('Fraction among singles')

    P_m_S_fig = figure('Color','White');
    plot(values,P_m_S_vec);
    legend({'No schooling','Primary','Lower middle','Upper middle',...
            'College'},'Location','best');
    legend('boxoff')
    xlabel(title_x_axis,'interpreter','latex')
    ylabel('Fraction among singles')

elseif agent_types == 2
    
    P_f_S_fig = figure('Color','White');
    plot(values,P_f_S_vec([10,30,50,70,90],:));
    legend({'p10','p30','p50','p70','p90'},'Location','best');
    legend('boxoff')
    xlabel(title_x_axis,'interpreter','latex')
    ylabel('Fraction among singles')

    P_m_S_fig = figure('Color','White');
    plot(values,P_m_S_vec([10,30,50,70,90],:));
    legend({'p10','p30','p50','p70','p90'},'Location','best');
    legend('boxoff')
    xlabel(title_x_axis,'interpreter','latex')
    ylabel('Fraction among singles')
    
end

match_prob_fig = figure('Color','White');
plot(values,[pi_f_vec; pi_m_vec]);
legend({'Females','Males'},'Location','northwest')
legend('boxoff')
xlabel(title_x_axis,'interpreter','latex')
ylabel('Matching probability')

if agent_types == 1
    
    MR_f_fig = figure('Color','White');
    plot(values,MR_f);
    legend({'No schooling','Primary','Lower middle','Upper middle',...
            'College'},'Location','best');
    legend('boxoff')
    xlabel(title_x_axis,'interpreter','latex')
    ylabel('Marriage rate per period')

    MR_m_fig = figure('Color','White');
    plot(values,MR_m);
    legend({'No schooling','Primary','Lower middle','Upper middle',...
            'College'},'Location','best');
    legend('boxoff')
    xlabel(title_x_axis,'interpreter','latex')
    ylabel('Marriage rate per period')
    
elseif agent_types == 2
    
    MR_f_fig = figure('Color','White');
    plot(values,MR_f([10,30,50,70,90],:));
    legend({'p10','p30','p50','p70','p90'},'Location','best');
    legend('boxoff')
    xlabel(title_x_axis,'interpreter','latex')
    ylabel('Marriage rate per period')

    P_m_S_fig = figure('Color','White');
    plot(values,MR_m([10,30,50,70,90],:));
    legend({'p10','p30','p50','p70','p90'},'Location','best');
    legend('boxoff')
    xlabel(title_x_axis,'interpreter','latex')
    ylabel('Marriage rate per period')
    
end

marr_rates_fig = figure('Color','White');
plot(values,[mr_f; mr_m]);
legend({'Females','Males'},'Location','northwest')
legend('boxoff')
xlabel(title_x_axis,'interpreter','latex')
ylabel('Marriage rate per period')

exp_single_fig = figure('Color','White');
plot(values,[exp_single_f; exp_single_m]);
legend({'Females','Males'},'Location','best')
legend('boxoff')
xlabel(title_x_axis,'interpreter','latex')
ylabel('Expected number of periods spent single')

sorting_fig = figure('Color','White');
plot(values,delta_vec);
legend('boxoff')
xlabel(title_x_axis,'interpreter','latex')
ylabel('\delta')

hypergamy_fig = figure('Color','White');
plot(values,[hypgam_f_vec; hypgam_m_vec])
legend({'Female', 'Male'},'Location','best')
legend('boxoff')
xlabel(title_x_axis,'interpreter','latex')
ylabel('Hypergamy measure')

sorting_2_fig = figure('Color','White');
plot(values,delta_2_vec);
legend('boxoff')
xlabel(title_x_axis,'interpreter','latex')
ylabel('Alternative \delta')

hypergamy_2_fig = figure('Color','White');
plot(values,[hypgam_f_2_vec; hypgam_m_2_vec])
legend({'Female', 'Male'},'Location','best')
legend('boxoff')
xlabel(title_x_axis,'interpreter','latex')
ylabel('Alternative hypergamy measure')

lf_fig = figure('Color','White');
plot(values,lf_bar);
legend('boxoff')
xlabel(title_x_axis,'interpreter','latex')
ylabel('Married female market hours as fraction of husbands')

ineq_fig = figure('Color','White');
plot(values,[ineq_inc; ineq_cons]);
legend({'Income','Consumption'},'Location','best')
legend('boxoff')
xlabel(title_x_axis,'interpreter','latex')
ylabel('Gini index')

if strcmp(save_graphs,'Yes')==1
    
    cd('Comparative statics/Graphs')
    
    saveas(VS_f_fig,strcat('VS_f_',graph_suffix),'epsc');
    saveas(VS_m_fig,strcat('VS_m_',graph_suffix),'epsc');
    saveas(match_prob_fig,strcat('match_prob_',graph_suffix),'epsc');
    saveas(MR_f_fig,strcat('MR_f_',graph_suffix),'epsc');
    saveas(MR_m_fig,strcat('MR_m_',graph_suffix),'epsc');
    saveas(marr_rates_fig,strcat('marr_rates_',graph_suffix),'epsc');
    saveas(sorting_fig,strcat('sorting_',graph_suffix),'epsc');
    saveas(hypergamy_fig,strcat('hypergamy_',graph_suffix),'epsc');
    saveas(sorting_2_fig,strcat('sorting_2_',graph_suffix),'epsc');
    saveas(hypergamy_2_fig,strcat('hypergamy_2_',graph_suffix),'epsc');
    saveas(exp_single_fig,strcat('exp_single_',graph_suffix),'epsc');
    saveas(lf_fig,strcat('lf_',graph_suffix),'epsc');
    saveas(ineq_fig,strcat('ineq_',graph_suffix),'epsc');
    
end



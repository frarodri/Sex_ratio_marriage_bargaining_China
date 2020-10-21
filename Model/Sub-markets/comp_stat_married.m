%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Comparative statics for the time allocation problem faced by married
% households


%%
% Preliminaries

% home: G:\My Drive\China project\Model\Sub-markets
% uni: G:\Mi unidad\China project\Model\Sub-markets

clear ; close all; clc
cd('G:\My Drive\China project\Model\Sub-markets')

%%
% Choose settings for the comparative statics exercise
  
promt = 'Primary parameter for comparative statics: ';
first = input(promt,'s');

promt2 = 'Minimum value: ';
min_first = input(promt2);

promt3 = 'Step: ';
step_first = input(promt3);

promt4 = 'Maximum value: ';
max_first = input(promt4);

promt5 = 'Title for x-axis in graphs: ';
title_x_axis = input(promt5,'s');

% promt5 = 'Perform comparative statics on a second parameter?';
% secondary_comp_stat = questdlg(promt5,'Secondary parameter');
% 
% if strcmp(secondary_comp_stat,'Yes')==1
%     
%     promt6 = 'Secondary parameter for comparative statics: ';
%     second = input(promt6,'s');
%     
%     promt7 = 'Vector of values for second: ';
%     vals2 = input(promt7,'s');
%     
% end

%%
% Baseline parameters

% Parameters of the utility function
param.sigma = 0.5;
param.sigma_c = 0.343;
param.sigma_l = 0.618;
param.sigma_g = 0.039;

% Parameters of the home production function
param.A_g = 1;
param.alpha_g = 0.95;

% Parameters of the housework aggregator
param.eta_f = 0.5;
param.eta = 0.33;

% Parameters that define wages
param.gender_wage_ratio = 0.75;
param.male_wage= 1;

% Parameter for the Pareto weight for the wife
param.pwf = 0.5;

% Parameter for price of home equipment
param.pe = 1;

%%
% Pre-allocate vectors

values_first = (min_first:step_first:max_first);
% values_second = str2num(vals2);

sz_first = size(values_first,2);
% sz_second = size(values_second);
sz_second = 1;

cf_vec = zeros(sz_first,sz_second);
cm_vec = zeros(sz_first,sz_second);
lf_vec = zeros(sz_first,sz_second);
lm_vec = zeros(sz_first,sz_second);
g_vec = zeros(sz_first,sz_second);
hf_vec = zeros(sz_first,sz_second);
hm_vec = zeros(sz_first,sz_second);
h_vec = zeros(sz_first,sz_second);
eq_vec = zeros(sz_first,sz_second);
nf_vec = zeros(sz_first,sz_second);
nm_vec = zeros(sz_first,sz_second);

%%
% Compute the solutions for each parameter value

for i = 1:sz_first
   
    param.(first) = values_first(i);
    wages_m = param.male_wage;
    wages_f = param.male_wage*param.gender_wage_ratio;
    pwf = param.pwf;
    
    [cf_vec(i),cm_vec(i),lf_vec(i),lm_vec(i),g_vec(i),hf_vec(i),...
        hm_vec(i),h_vec(i),eq_vec(i),nf_vec(i),nm_vec(i)]=...
        per_sol_married(param,wages_f,wages_m,pwf,'quietly','true');
    
end

%%
% Graphs

cf_fig = figure('Color','White');
plot(values_first,cf_vec);
xlabel(title_x_axis,'interpreter','latex')
ylabel("Wife's consumption")

cm_fig = figure('Color','White');
plot(values_first,cm_vec);
xlabel(title_x_axis,'interpreter','latex')
ylabel("Husband's consumption")

lf_fig = figure('Color','White');
plot(values_first,lf_vec);
xlabel(title_x_axis,'interpreter','latex')
ylabel("Wife's paid work time")

lm_fig = figure('Color','White');
plot(values_first,lm_vec);
xlabel(title_x_axis,'interpreter','latex')
ylabel("Husband's paid work time")

g_fig = figure('Color','White');
plot(values_first,g_vec);
xlabel(title_x_axis,'interpreter','latex')
ylabel("Home-produced good consumption")

hf_fig = figure('Color','White');
plot(values_first,hf_vec);
xlabel(title_x_axis,'interpreter','latex')
ylabel("Wife's housework time")

hm_fig = figure('Color','White');
plot(values_first,hm_vec);
xlabel(title_x_axis,'interpreter','latex')
ylabel("Husband's housework time")

h_fig = figure('Color','White');
plot(values_first,h_vec);
xlabel(title_x_axis,'interpreter','latex')
ylabel("Aggregate housework time")

eq_fig = figure('Color','White');
plot(values_first,eq_vec);
xlabel(title_x_axis,'interpreter','latex')
ylabel("Equipment usage")

nf_fig = figure('Color','White');
plot(values_first,nf_vec);
xlabel(title_x_axis,'interpreter','latex')
ylabel("Wife's leisure time")

nm_fig = figure('Color','White');
plot(values_first,nm_vec);
xlabel(title_x_axis,'interpreter','latex')
ylabel("Husband's leisure time")
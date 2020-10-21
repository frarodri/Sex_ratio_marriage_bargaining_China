%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Comparative statics for the time allocation problem faced by married
% households


%%
% Preliminaries

% home: G:\My Drive\China project\Model\Last version
% uni: G:\Mi unidad\China project\Model\Sub-markets

clear ; close all; clc
cd('G:\My Drive\China project\Model\Last version')

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

promt6 = 'Perform comparative statics on a second parameter?';
secondary_comp_stat = questdlg(promt6,'Secondary parameter');
 
if strcmp(secondary_comp_stat,'Yes')==1
    
    promt7 = 'Secondary parameter for comparative statics: ';
    second = input(promt7,'s');
    
    promt8 = 'Vector of values for second: ';
    vals2 = input(promt8,'s');
    
    promt9 = 'Variable name in legend: ';
    var_legend = input(promt9,'s');
    
end

%%
% Baseline parameters

% Parameters of the utility function
param.sigma = 1.5;
param.sigma_c = 0.3533;
param.sigma_l = 0.6055;
param.sigma_g = 1-param.sigma_c-param.sigma_l;

% Parameters of the home production function
param.A_g = 1;
param.alpha_g = 0.95;

% Parameters of the housework aggregator
param.eta_f = 0.5801;
param.eta = 0.33;

% Parameters that define wages
param.gender_wage_ratio = 0.83;
param.male_wage= 1;

% Parameter for price of home equipment
param.pe = 1.5;

% Parameter for the Pareto weight for the wife
param.pwf = 0.5;

%%
% Pre-allocate vectors

values_first = (min_first:step_first:max_first);
sz_first = size(values_first,2);

if strcmp(secondary_comp_stat,'Yes')==1

    values_second = str2num(vals2);
    sz_second = size(values_second,2);
    
    legend_text = string(var_legend)+"="+values_second;
    
else
   
    sz_second = 1;
    
end

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
% Compute the solutions for each set of parameter values

if strcmp(secondary_comp_stat,'Yes')==1
    
    for i = 1:sz_first
        
        for j = 1:sz_second
            
            param.(first) = values_first(i);
            param.(second) = values_second(j);
            
            wages_m = param.male_wage;
            wages_f = param.male_wage*param.gender_wage_ratio;
            pwf = param.pwf;
            
            [cf_vec(i,j),cm_vec(i,j),lf_vec(i,j),lm_vec(i,j),g_vec(i,j),...
                hf_vec(i,j),hm_vec(i,j),h_vec(i,j),eq_vec(i,j),...
                nf_vec(i,j),nm_vec(i,j)]=per_sol_married(param,wages_f,...
                wages_m,pwf,'quietly','true');
            
        end
        
    end
    
    % Graphs

    cf_fig = figure('Color','White');
    plot(values_first,cf_vec);
    legend(legend_text,'Interpreter','latex','Location','best')
    legend('boxoff')
    xlabel(title_x_axis,'interpreter','latex')
    ylabel("Wife's consumption")

    cm_fig = figure('Color','White');
    plot(values_first,cm_vec);
    legend(legend_text,'Interpreter','latex','Location','best')
    legend('boxoff')
    xlabel(title_x_axis,'interpreter','latex')
    ylabel("Husband's consumption")

    lf_fig = figure('Color','White');
    plot(values_first,lf_vec);
    legend(legend_text,'Interpreter','latex','Location','best')
    legend('boxoff')
    xlabel(title_x_axis,'interpreter','latex')
    ylabel("Wife's leisure time")

    lm_fig = figure('Color','White');
    plot(values_first,lm_vec);
    legend(legend_text,'Interpreter','latex','Location','best')
    legend('boxoff')
    xlabel(title_x_axis,'interpreter','latex')
    ylabel("Husband's leisure time")

    g_fig = figure('Color','White');
    plot(values_first,g_vec);
    legend(legend_text,'Interpreter','latex','Location','best')
    legend('boxoff')
    xlabel(title_x_axis,'interpreter','latex')
    ylabel("Home-produced good consumption")

    hf_fig = figure('Color','White');
    plot(values_first,hf_vec);
    legend(legend_text,'Interpreter','latex','Location','best')
    legend('boxoff')
    xlabel(title_x_axis,'interpreter','latex')
    ylabel("Wife's housework time")

    hm_fig = figure('Color','White');
    plot(values_first,hm_vec);
    legend(legend_text,'Interpreter','latex','Location','best')
    legend('boxoff')
    xlabel(title_x_axis,'interpreter','latex')
    ylabel("Husband's housework time")

    h_fig = figure('Color','White');
    plot(values_first,h_vec);
    legend(legend_text,'Interpreter','latex','Location','best')
    legend('boxoff')
    xlabel(title_x_axis,'interpreter','latex')
    ylabel("Aggregate housework time")

    eq_fig = figure('Color','White');
    plot(values_first,eq_vec);
    legend(legend_text,'Interpreter','latex','Location','best')
    legend('boxoff')
    xlabel(title_x_axis,'interpreter','latex')
    ylabel("Equipment usage")

    nf_fig = figure('Color','White');
    plot(values_first,nf_vec);
    legend(legend_text,'Interpreter','latex','Location','best')
    legend('boxoff')
    xlabel(title_x_axis,'interpreter','latex')
    ylabel("Wife's paid work time")

    nm_fig = figure('Color','White');
    plot(values_first,nm_vec);
    legend(legend_text,'Interpreter','latex','Location','best')
    legend('boxoff')
    xlabel(title_x_axis,'interpreter','latex')
    ylabel("Husband's paid work time")
    
else

    for i = 1:sz_first

        param.(first) = values_first(i);
        wages_m = param.male_wage;
        wages_f = param.male_wage*param.gender_wage_ratio;
        pwf = param.pwf;

        [cf_vec(i),cm_vec(i),lf_vec(i),lm_vec(i),g_vec(i),hf_vec(i),...
            hm_vec(i),h_vec(i),eq_vec(i),nf_vec(i),nm_vec(i)]=...
            per_sol_married(param,wages_f,wages_m,pwf,'quietly','true');

    end

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
    ylabel("Wife's leisure time")

    lm_fig = figure('Color','White');
    plot(values_first,lm_vec);
    xlabel(title_x_axis,'interpreter','latex')
    ylabel("Husband's leisure time")

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
    ylabel("Wife's paid work time")

    nm_fig = figure('Color','White');
    plot(values_first,nm_vec);
    xlabel(title_x_axis,'interpreter','latex')
    ylabel("Husband's paid work time")
    
end
%{
This code calibrates the model to the data, performs the quantitative
experiments and organizes the results 
%} 

%% Preliminaries
clear ; close all; clc
main_dir = 'G:\My Drive\China project\Model\Last version';
data_dir = 'G:\My Drive\China project\Data\JMK';
results_dir = 'G:\My Drive\China project\Paper\JMK';

filename_data = 'Data_for_model.xlsx';
filename_results = 'Results.xlsx';

%% Choose sigma and whether to calibrate or not

promt = 'Value for sigma: ';
value_sigma = input(promt);

promt2 = 'Perform calibration?';
calibrate = questdlg(promt2,'Calibrate?');

%% Read the data

cd(data_dir);

filename = 'Data_for_model.xlsx';

% Sex ratio
T = readtable(filename_data,'Sheet','sexratio');
data.sexratio = table2struct(T);

% Time allocation data
T = readtable(filename_data,'Sheet','ta1990');
data.ta1990 = table2struct(T);

T = readtable(filename_data,'Sheet','ta2010');
data.ta2010 = table2struct(T);

% Skill ditributions
T =  readtable(filename_data,'Sheet','skilldist');
data.skilldist = table2struct(T,'ToScalar',true);

% Marital sorting data
data.MS.abs1990 = readmatrix(filename_data,...
    'Sheet','Marital Sorting 1990',...
    'Range','B2');

data.MS.abs2010 = readmatrix(filename_data,...
    'Sheet','Marital Sorting 2010',...
    'Range','B2');

% Wages
T = readtable(filename_data,'Sheet','wages1990');
data.wages1990 = table2struct(T);

T = readtable(filename_data,'Sheet','wages2010');
data.wages2010 = table2struct(T);

%% Perform the appropriate transformations of the data

cd(main_dir);

% Marital sorting
% Obtain the contingency matrices
data.MS.CM1990 = data.MS.abs1990/sum(data.MS.abs1990,'all');
data.MS.CM2010 = data.MS.abs2010/sum(data.MS.abs2010,'all');
% Obtain measures of marital sorting
[data.MS.delta1990, data.MS.hypf1990, data.MS.hypm1990] = ...
    msmeasures(data.MS.CM1990);
[data.MS.delta2010, data.MS.hypf2010, data.MS.hypm2010] = ...
    msmeasures(data.MS.CM2010);

% Wages
wages.m1990 = [data.wages1990.lowskill;...
    data.wages1990.midskill;data.wages1990.highskill]/...
    data.wages1990.lowskill;

wages.m2010 = [data.wages2010.lowskill;...
    data.wages2010.midskill;data.wages2010.highskill]/...
    data.wages1990.lowskill;

wages.phi1990 = data.wages1990.female/data.wages1990.male;
wages.phi2010 = data.wages2010.female/data.wages2010.male;

wages.f1990 = wages.m1990*wages.phi1990;
wages.f2010 = wages.m2010*wages.phi2010;

wages.skillpremium1990 = [data.wages1990.lowskill;...
    data.wages1990.midskill;data.wages1990.highskill]/...
    data.wages1990.lowskill;

wages.skillpremium2010 = [data.wages2010.lowskill;...
    data.wages2010.midskill;data.wages2010.highskill]/...
    data.wages2010.lowskill;

wages.growth = data.wages2010.all/data.wages1990.all;

% Export data on skill distributions and wages

cd(results_dir);

Sheet = "skills_distribution";

skills_distribution = struct2array(data.skilldist);

writematrix(skills_distribution,filename_results,'Sheet',Sheet);

Sheet = "wages";

wages_export = [wages.m1990,wages.f1990,wages.m2010,wages.f2010];

writematrix(wages_export,filename_results,'Sheet',Sheet);

Sheet = "sex ratio";

sex_ratio = struct2table(data.sexratio);

writetable(sex_ratio,filename_results,'Sheet',Sheet);

%% Parameters set externally

% Parameters that remain constant
param.beta = 0.96; % Standard
param.delta = 1/49; % Life expectancy of 69 years according to UN 
                    % Population Division in 1990
param.rho = 1/15; % Expected number of periods in marriage market
param.sigma = value_sigma; % Attanasio et. al (2008)
param.alpha_g = 1-0.05; % U.S. number in Knowles (2014)
param.eta = 0.33; % Knowles (2014)

ext_param = struct2table(param);

cd(results_dir);

Sheet = "external_parameters, sigma=" + param.sigma;

writetable(ext_param,filename_results,'Sheet',Sheet);

param.alpha_x = 1/2; % Makes sense for a matching function
param.A_x = 1;
param.growth_A_g = 0.0958; % Growth rate GDP per capita 1990-2010

%% Compute the remaining parameter of the housework time aggregator

housewrkrat1990 = data.ta1990.hfmarried/data.ta1990.hmmarried;
housewrkrat2010 = data.ta2010.hfmarried/data.ta2010.hmmarried;

x = (housewrkrat1990^param.eta)*wages.phi1990;

param.eta_f = x/(1+x);

% Check how well this adjusts to data
eta_implied = (log(1/wages.phi2010)-log(1/wages.phi1990))/...
    (log(housewrkrat2010)-log(housewrkrat1990));


%% Compute the weights in the utility function for single people

cd(main_dir);

% Set the price of home equipment and efficiency of home production to 1990
% values
param.pe = 1.5;
param.A_g = 1;

% Define the inputs for the minimization routine
x0 = 1/3*[1 1 1];
A = [];
b = []; 
Aeq = [1 1 1];
beq = 1;
lb = [0 0 0];
ub = [];

options = optimoptions('fmincon','Display','off');

% Single women
% Define the targets 
targ.h = data.ta1990.hfsingle;
targ.n = data.ta1990.nfsingle;
targ.l = data.ta1990.lfsingle;

% Set the wage
wage = sum(data.skilldist.Pf1990.*wages.f1990);

% Minimize the loss function
uw_singlefemale = fmincon(@(x)timealloclossfn_s(param,targ,wage,x),x0,A,...
                    b,Aeq,beq,lb,ub,[],options);

% Single men
% Define the targets
targ.h = data.ta1990.hmsingle;
targ.n = data.ta1990.nmsingle;
targ.l = data.ta1990.lmsingle;

% Set the wage
wage = sum(data.skilldist.Pm1990.*wages.m1990);

% Minimize the loss function
uw_singlemale = fmincon(@(x)timealloclossfn_s(param,targ,wage,x),x0,A,b,...
                                                Aeq,beq,lb,ub,[],options);
                                            
% Export the set of parameters that are chosen to match data before solving
% the model

data_match_param = [param.eta_f;uw_singlefemale';uw_singlemale'];

cd(results_dir);

Sheet = "data_match_param, sigma=" + param.sigma;

writematrix(data_match_param,filename_results,'Sheet',Sheet);
                                            
%% Compute the period utility flow for singles

cd(main_dir);

% Females
wage = wages.f1990;
param.sigma_c = uw_singlefemale(1);
param.sigma_l = uw_singlefemale(2);
param.sigma_g = uw_singlefemale(3);
USf1990 = uflow_singles(param,wage);

% Males
wage = wages.m1990;
param.sigma_c = uw_singlemale(1);
param.sigma_l = uw_singlemale(2);
param.sigma_g = uw_singlemale(3);
USm1990 = uflow_singles(param,wage);

%%
if strcmp(calibrate,'Yes')==1
    
    %% Calibration
    % Here we need to compute the utility weights for married couples, the 
    % cost of being single for females and means of the match quality draws 
    % that minimize the distance between the data targets and the model
    % counterparts. 

    cd(main_dir);

    % Set the parameters, wages and distributions to 1990 values
    param.pe = 1.5;
    param.A_g = 1;

    param.theta_0 = data.sexratio.sexratio1990; 

    param.Pf = data.skilldist.Pf1990;
    param.Pm = data.skilldist.Pm1990;

    wages.m = wages.m1990;
    wages.f = wages.f1990;

    targ.lf = data.ta1990.lfmarried;
    targ.lm = data.ta1990.lmmarried;
    targ.hf = data.ta1990.hfmarried;
    targ.hm = data.ta1990.hmmarried;
    targ.nf = data.ta1990.nfmarried;
    targ.nm = data.ta1990.nmmarried;

    targ.MS = data.MS.CM1990; 

    min_sigma_c = 0.35;
    max_sigma_c = 0.425;

    min_sigma_l = 0.50;
    max_sigma_l = 0.625;

    min_psi = -1;
    max_psi = 0;

    current_err = 1;
    tol = 10^(-4);
    it = 1;

    param.sigma_c = (max_sigma_c+min_sigma_c)/2;
    param.sigma_l = (max_sigma_l+min_sigma_l)/2;
    param.sigma_g = 1-param.sigma_c-param.sigma_l;

    param.psi = (max_psi+min_psi)/2;

    param.MU = zeros(size(param.Pf,1));

    t = cputime;

    while current_err>tol

        [lf,lm,hf,hm,nf,nm,mus] = ...
            calib_iteration(param,wages,targ,USf1990,USm1990);

        if nf+nm>targ.nf+targ.nm

            max_sigma_c = param.sigma_c;

        else

            min_sigma_c = param.sigma_c;

        end

        if lf+lm>targ.lf+targ.lm

            max_sigma_l = param.sigma_l;

        else

            min_sigma_l = param.sigma_l;

        end

        if lf/lm>targ.lf/targ.lm

            max_psi = param.psi;

        else

            min_psi = param.psi;

        end

        param.sigma_c = (max_sigma_c+min_sigma_c)/2;
        param.sigma_l = (max_sigma_l+min_sigma_l)/2;
        param.sigma_g = 1-param.sigma_c-param.sigma_l;

        param.psi = (max_psi+min_psi)/2;

        current_err = max(max_sigma_c-min_sigma_c,...
            max(max_sigma_l-min_sigma_l,max_psi-min_psi));

    end

    e = cputime-t;

    message_elapsed_time = 'Elapsed time for calibration: %g seconds \n';
    fprintf(message_elapsed_time,e);
    
    married_sigma_c = param.sigma_c;
    married_sigma_l = param.sigma_l;
    married_sigma_g = param.sigma_g;
    psi = param.psi;
    
    mus_mat = vec2mat(mus,size(param.Pf,1));
    param.MU = mus_mat;

    % Export the calibration results

    cd(results_dir);

    calib_res_ufn.Parameter = {'$\sigma_c$';'$\sigma_l$';'$\sigma_g$';...
        '$\psi_f$'};

    calib_res_ufn.Value = [param.sigma_c ; param.sigma_l; param.sigma_g;...
        param.psi];

    calibration_results_ufn = struct2table(calib_res_ufn);

    calib_res_ufn.Target = [married_sigma_c ; married_sigma_l; ...
                            married_sigma_g ; psi];

    Sheet = "calib_results_u, sigma=" + param.sigma;

    writetable(calibration_results_ufn,filename_results,'Sheet',Sheet);

    calib_res_ms.Woman_man = {'Low';'Medium';'High'};
    calib_res_ms.Low = mus_mat(:,1);
    calib_res_ms.Medium = mus_mat(:,2);
    calib_res_ms.High = mus_mat(:,3);

    calibration_results_ms = struct2table(calib_res_ms);

    Sheet = "calib_results_ms, sigma=" + param.sigma;

    writetable(calibration_results_ms,filename_results,'Sheet',Sheet);
    %%
else 

    %% Set the calibrated parameters (when we dont want to do calibration)

    cd(results_dir);
    Sheet = "Calib_results_u, sigma=" + param.sigma;

    married_utility_parameters = ...
        readmatrix(filename_results,...
        'Sheet',Sheet,...
        'Range','B2');

    married_sigma_c = married_utility_parameters(1);
    married_sigma_l = married_utility_parameters(2);
    married_sigma_g = married_utility_parameters(3);
    psi_f = married_utility_parameters(4);

    param.sigma_c = married_sigma_c;
    param.sigma_l = married_sigma_l;
    param.sigma_g = married_sigma_g;
    param.psi = psi_f;

    Sheet = "Calib_results_ms, sigma=" + param.sigma;

    param.MU = readmatrix(filename_results,...
        'Sheet',Sheet,...
        'Range','B2');
    %%
end

%% Compute and organize the results for the baseline year (1990)

cd(main_dir);

% Set the parameters, wages and distributions to 1990 values
param.pe = 1.5;
param.A_g = 1;

param.theta_0 = data.sexratio.sexratio1990;

param.Pf = data.skilldist.Pf1990;
param.Pm = data.skilldist.Pm1990;

wages.m = wages.m1990;
wages.f = wages.f1990;

% Compute equilibrium marriage objects
[theta_S_1990,Q_1990,PWf_1990] = sseq(param,wages,USf1990,USm1990,...
    'bargaining','Egalitarian');

MS_1990 = marital_sorting(param,theta_S_1990,Q_1990);
model.MS.CM1990 = MS_1990;

avgPWf1990 = sum(sum(MS_1990.*PWf_1990));

[Cf_1990,Cm_1990,Lf_1990,Lm_1990,G_1990,Hf_1990,Hm_1990,H_1990,EQ_1990,...
    Nf_1990,Nm_1990] = per_sol_married(param,wages.f,wages.m,PWf_1990);
            
model.ta1990.hfmarried = sum(sum(MS_1990.*Hf_1990));
model.ta1990.nfmarried = sum(sum(MS_1990.*Nf_1990));
model.ta1990.lfmarried = sum(sum(MS_1990.*Lf_1990));
model.ta1990.hmmarried = sum(sum(MS_1990.*Hm_1990));
model.ta1990.nmmarried = sum(sum(MS_1990.*Nm_1990));
model.ta1990.lmmarried = sum(sum(MS_1990.*Lm_1990));

% Single female time allocation
param.sigma_c = uw_singlefemale(1);
param.sigma_l = uw_singlefemale(2);
param.sigma_g = uw_singlefemale(3);
[~,lfs,~,hfs,~,nfs] = per_sol_singles(param,wages.f);

model.ta1990.hfsingle = sum(hfs.*param.Pf);
model.ta1990.nfsingle = sum(nfs.*param.Pf);
model.ta1990.lfsingle = sum(lfs.*param.Pf);

% Single men time allocation
param.sigma_c = uw_singlemale(1);
param.sigma_l = uw_singlemale(2);
param.sigma_g = uw_singlemale(3);

[~,lms,~,hms,~,nms] = per_sol_singles(param,wages.m);

model.ta1990.hmsingle = sum(hms.*param.Pm);
model.ta1990.nmsingle = sum(nms.*param.Pm);
model.ta1990.lmsingle = sum(lms.*param.Pm);

% Export results
cd(results_dir);
% Table comparing time allocation in the model and in the data

ta_results_1990.Statistic = {'Married women housework';...
    'Married women paid work';'Married women leisure';...
    'Married men housework';'Married men paid work';...
    'Married men leisure';'Single women housework';...
    'Single women paid work';'Single women leisure';...
    'Single men housework';'Single men paid work';...
    'Single men leisure'};

ta_results_1990.Data = struct2array(data.ta1990)'*118;
ta_results_1990.Model = struct2array(model.ta1990)'*118;

time_allocation_results_1990 = struct2table(ta_results_1990);

Sheet = "ta_res_1990, sigma=" + param.sigma;

writetable(time_allocation_results_1990,filename_results,'Sheet',Sheet);

% Table comparing the marital sorting in the model and in the data

ms_results_1990.Wife_husband = {'Low skill';'Medium skill';'High skill';...
    'Low skill';'Medium skill';'High skill'};
ms_results_1990.Low_skill = [data.MS.CM1990(:,1);model.MS.CM1990(:,1)];
ms_results_1990.Medium_skill = [data.MS.CM1990(:,2);model.MS.CM1990(:,2)];
ms_results_1990.High_skill = [data.MS.CM1990(:,3);model.MS.CM1990(:,3)];

marital_sorting_results_1990 = struct2table(ms_results_1990);

Sheet = "ms_res_1990, sigma=" + param.sigma;

writetable(marital_sorting_results_1990,filename_results,'Sheet',Sheet);


%% Compute and organize the results for 2010

cd(main_dir);

% Set the parameters, wages and distributions to 2010 values
param.pe = 0.5;
param.A_g = (1+param.growth_A_g)^20;

param.theta_0 = data.sexratio.sexratio2010;

param.Pf = data.skilldist.Pf2010;
param.Pm = data.skilldist.Pm2010;

wages.m = wages.m2010;
wages.f = wages.f2010;

param.sigma_c = uw_singlefemale(1);
param.sigma_l = uw_singlefemale(2);
param.sigma_g = uw_singlefemale(3);
USf2010 = uflow_singles(param,wages.f);

param.sigma_c = uw_singlemale(1);
param.sigma_l = uw_singlemale(2);
param.sigma_g = uw_singlemale(3);
USm2010 = uflow_singles(param,wages.m);

param.sigma_c = married_sigma_c;
param.sigma_l = married_sigma_l;
param.sigma_g = 1-param.sigma_c-param.sigma_l;

% Compute equilibrium marriage objects
[theta_S_2010,Q_2010,PWf_2010] = sseq(param,wages,USf2010,USm2010,...
    'bargaining','Egalitarian');

MS_2010 = marital_sorting(param,theta_S_2010,Q_2010);
model.MS.CM2010 = MS_2010;

avgPWf2010 = sum(sum(MS_2010.*PWf_2010));

[Cf_2010,Cm_2010,Lf_2010,Lm_2010,G_2010,Hf_2010,Hm_2010,H_2010,EQ_2010,...
    Nf_2010,Nm_2010] = per_sol_married(param,wages.f,wages.m,PWf_2010);
            
model.ta2010.hfmarried = sum(sum(MS_2010.*Hf_2010));
model.ta2010.nfmarried = sum(sum(MS_2010.*Nf_2010));
model.ta2010.lfmarried = sum(sum(MS_2010.*Lf_2010));
model.ta2010.hmmarried = sum(sum(MS_2010.*Hm_2010));
model.ta2010.nmmarried = sum(sum(MS_2010.*Nm_2010));
model.ta2010.lmmarried = sum(sum(MS_2010.*Lm_2010));

% Single female time allocation
param.sigma_c = uw_singlefemale(1);
param.sigma_l = uw_singlefemale(2);
param.sigma_g = uw_singlefemale(3);
[~,lfs,~,hfs,~,nfs] = per_sol_singles(param,wages.f);

model.ta2010.hfsingle = sum(hfs.*param.Pf);
model.ta2010.nfsingle = sum(nfs.*param.Pf);
model.ta2010.lfsingle = sum(lfs.*param.Pf);

% Single men time allocation
param.sigma_c = uw_singlemale(1);
param.sigma_l = uw_singlemale(2);
param.sigma_g = uw_singlemale(3);

[~,lms,~,hms,~,nms] = per_sol_singles(param,wages.m);

model.ta2010.hmsingle = sum(hms.*param.Pm);
model.ta2010.nmsingle = sum(nms.*param.Pm);
model.ta2010.lmsingle = sum(lms.*param.Pm);

% Export results
cd(results_dir);
% Table comparing time allocation in the model and in the data

ta_results_2010.Statistic = {'Married women housework';...
    'Married women paid work';'Married women leisure';...
    'Married men housework';'Married men paid work';...
    'Married men leisure';'Single women housework';...
    'Single women paid work';'Single women leisure';...
    'Single men housework';'Single men paid work';...
    'Single men leisure'};

ta_results_2010.Data = struct2array(data.ta2010)'*118;
ta_results_2010.Model = struct2array(model.ta2010)'*118;

time_allocation_results_2010 = struct2table(ta_results_2010);

Sheet = "ta_res_2010, sigma=" + param.sigma;

writetable(time_allocation_results_2010,filename_results,'Sheet',Sheet);

% Table comparing the marital sorting in the model and in the data

ms_results_2010.Wife_husband = {'Low skill';'Medium skill';'High skill';...
    'Low skill';'Medium skill';'High skill'};
ms_results_2010.Low_skill = [data.MS.CM2010(:,1);model.MS.CM2010(:,1)];
ms_results_2010.Medium_skill = [data.MS.CM2010(:,2);model.MS.CM2010(:,2)];
ms_results_2010.High_skill = [data.MS.CM2010(:,3);model.MS.CM2010(:,3)];

marital_sorting_results_2010 = struct2table(ms_results_2010);

Sheet = "ms_res_2010, sigma=" + param.sigma;

writetable(marital_sorting_results_2010,filename_results,'Sheet',Sheet);

%% Decomposition (changing parameters one at a time)

% Backward decomposition

%% First stage: the sex ratio goes back to 1990 level

cd(main_dir);

param.theta_0 = data.sexratio.sexratio1990;

param.sigma_c = uw_singlefemale(1);
param.sigma_l = uw_singlefemale(2);
param.sigma_g = uw_singlefemale(3);
USfd1 = uflow_singles(param,wages.f);

param.sigma_c = uw_singlemale(1);
param.sigma_l = uw_singlemale(2);
param.sigma_g = uw_singlemale(3);
USmd1 = uflow_singles(param,wages.m);

param.sigma_c = married_sigma_c;
param.sigma_l = married_sigma_l;
param.sigma_g = 1-param.sigma_c-param.sigma_l;

% Compute equilibrium marriage objects
[theta_S_d1,Q_d1,PWf_d1] = sseq(param,wages,USfd1,USmd1,...
    'bargaining','Egalitarian');

MS_d1 = marital_sorting(param,theta_S_d1,Q_d1);
decomp.MS.CMd1 = MS_d1;

avgPWfd1 = sum(sum(MS_d1.*PWf_d1));

[Cf_d1,Cm_d1,Lf_d1,Lm_d1,G_d1,Hf_d1,Hm_d1,H_d1,EQ_d1,...
    Nf_d1,Nm_d1] = per_sol_married(param,wages.f,wages.m,PWf_d1);
            
decomp.tad1.hfmarried = sum(sum(MS_d1.*Hf_d1));
decomp.tad1.nfmarried = sum(sum(MS_d1.*Nf_d1));
decomp.tad1.lfmarried = sum(sum(MS_d1.*Lf_d1));
decomp.tad1.hmmarried = sum(sum(MS_d1.*Hm_d1));
decomp.tad1.nmmarried = sum(sum(MS_d1.*Nm_d1));
decomp.tad1.lmmarried = sum(sum(MS_d1.*Lm_d1));

% Single female time allocation
param.sigma_c = uw_singlefemale(1);
param.sigma_l = uw_singlefemale(2);
param.sigma_g = uw_singlefemale(3);
[~,lfs,~,hfs,~,nfs] = per_sol_singles(param,wages.f);

decomp.tad1.hfsingle = sum(hfs.*param.Pf);
decomp.tad1.nfsingle = sum(nfs.*param.Pf);
decomp.tad1.lfsingle = sum(lfs.*param.Pf);

% Single men time allocation
param.sigma_c = uw_singlemale(1);
param.sigma_l = uw_singlemale(2);
param.sigma_g = uw_singlemale(3);

[~,lms,~,hms,~,nms] = per_sol_singles(param,wages.m);

decomp.tad1.hmsingle = sum(hms.*param.Pm);
decomp.tad1.nmsingle = sum(nms.*param.Pm);
decomp.tad1.lmsingle = sum(lms.*param.Pm);

%% Second stage: the distributions of skill go back to 1990 level

param.Pf = data.skilldist.Pf1990;
param.Pm = data.skilldist.Pm1990;

param.sigma_c = uw_singlefemale(1);
param.sigma_l = uw_singlefemale(2);
param.sigma_g = uw_singlefemale(3);
USfd2 = uflow_singles(param,wages.f);

param.sigma_c = uw_singlemale(1);
param.sigma_l = uw_singlemale(2);
param.sigma_g = uw_singlemale(3);
USmd2 = uflow_singles(param,wages.m);

param.sigma_c = married_sigma_c;
param.sigma_l = married_sigma_l;
param.sigma_g = 1-param.sigma_c-param.sigma_l;

% Compute equilibrium marriage objects
[theta_S_d2,Q_d2,PWf_d2] = sseq(param,wages,USfd2,USmd2,...
    'bargaining','Egalitarian');

MS_d2 = marital_sorting(param,theta_S_d2,Q_d2);
decomp.MS.CMd2 = MS_d2;

avgPWfd2 = sum(sum(MS_d2.*PWf_d2));

[Cf_d2,Cm_d2,Lf_d2,Lm_d2,G_d2,Hf_d2,Hm_d2,H_d2,EQ_d2,...
    Nf_d2,Nm_d2] = per_sol_married(param,wages.f,wages.m,PWf_d2);
            
decomp.tad2.hfmarried = sum(sum(MS_d2.*Hf_d2));
decomp.tad2.nfmarried = sum(sum(MS_d2.*Nf_d2));
decomp.tad2.lfmarried = sum(sum(MS_d2.*Lf_d2));
decomp.tad2.hmmarried = sum(sum(MS_d2.*Hm_d2));
decomp.tad2.nmmarried = sum(sum(MS_d2.*Nm_d2));
decomp.tad2.lmmarried = sum(sum(MS_d2.*Lm_d2));

% Single female time allocation
param.sigma_c = uw_singlefemale(1);
param.sigma_l = uw_singlefemale(2);
param.sigma_g = uw_singlefemale(3);
[~,lfs,~,hfs,~,nfs] = per_sol_singles(param,wages.f);

decomp.tad2.hfsingle = sum(hfs.*param.Pf);
decomp.tad2.nfsingle = sum(nfs.*param.Pf);
decomp.tad2.lfsingle = sum(lfs.*param.Pf);

% Single men time allocation
param.sigma_c = uw_singlemale(1);
param.sigma_l = uw_singlemale(2);
param.sigma_g = uw_singlemale(3);

[~,lms,~,hms,~,nms] = per_sol_singles(param,wages.m);

decomp.tad2.hmsingle = sum(hms.*param.Pm);
decomp.tad2.nmsingle = sum(nms.*param.Pm);
decomp.tad2.lmsingle = sum(lms.*param.Pm);

%% Third stage: the wage structure go back to 1990 level

wages.m = wages.m1990;
wages.f = wages.f1990;

param.sigma_c = uw_singlefemale(1);
param.sigma_l = uw_singlefemale(2);
param.sigma_g = uw_singlefemale(3);
USfd3 = uflow_singles(param,wages.f);

param.sigma_c = uw_singlemale(1);
param.sigma_l = uw_singlemale(2);
param.sigma_g = uw_singlemale(3);
USmd3 = uflow_singles(param,wages.m);

param.sigma_c = married_sigma_c;
param.sigma_l = married_sigma_l;
param.sigma_g = 1-param.sigma_c-param.sigma_l;

% Compute equilibrium marriage objects
[theta_S_d3,Q_d3,PWf_d3] = sseq(param,wages,USfd3,USmd3,...
    'bargaining','Egalitarian');

MS_d3 = marital_sorting(param,theta_S_d3,Q_d3);
decomp.MS.CMd3 = MS_d3;

avgPWfd3 = sum(sum(MS_d3.*PWf_d3));

[Cf_d3,Cm_d3,Lf_d3,Lm_d3,G_d3,Hf_d3,Hm_d3,H_d3,EQ_d3,...
    Nf_d3,Nm_d3] = per_sol_married(param,wages.f,wages.m,PWf_d3);
            
decomp.tad3.hfmarried = sum(sum(MS_d3.*Hf_d3));
decomp.tad3.nfmarried = sum(sum(MS_d3.*Nf_d3));
decomp.tad3.lfmarried = sum(sum(MS_d3.*Lf_d3));
decomp.tad3.hmmarried = sum(sum(MS_d3.*Hm_d3));
decomp.tad3.nmmarried = sum(sum(MS_d3.*Nm_d3));
decomp.tad3.lmmarried = sum(sum(MS_d3.*Lm_d3));

% Single female time allocation
param.sigma_c = uw_singlefemale(1);
param.sigma_l = uw_singlefemale(2);
param.sigma_g = uw_singlefemale(3);
[~,lfs,~,hfs,~,nfs] = per_sol_singles(param,wages.f);

decomp.tad3.hfsingle = sum(hfs.*param.Pf);
decomp.tad3.nfsingle = sum(nfs.*param.Pf);
decomp.tad3.lfsingle = sum(lfs.*param.Pf);

% Single men time allocation
param.sigma_c = uw_singlemale(1);
param.sigma_l = uw_singlemale(2);
param.sigma_g = uw_singlemale(3);

[~,lms,~,hms,~,nms] = per_sol_singles(param,wages.m);

decomp.tad3.hmsingle = sum(hms.*param.Pm);
decomp.tad3.nmsingle = sum(nms.*param.Pm);
decomp.tad3.lmsingle = sum(lms.*param.Pm);

%% Fourth stage: pe and Ag go back to 1990 level

param.pe = 1.5;
param.A_g = 1;

param.sigma_c = uw_singlefemale(1);
param.sigma_l = uw_singlefemale(2);
param.sigma_g = uw_singlefemale(3);
USfd4 = uflow_singles(param,wages.f);

param.sigma_c = uw_singlemale(1);
param.sigma_l = uw_singlemale(2);
param.sigma_g = uw_singlemale(3);
USmd4 = uflow_singles(param,wages.m);

param.sigma_c = married_sigma_c;
param.sigma_l = married_sigma_l;
param.sigma_g = 1-param.sigma_c-param.sigma_l;

% Compute equilibrium marriage objects
[theta_S_d4,Q_d4,PWf_d4] = sseq(param,wages,USfd4,USmd4,...
    'bargaining','Egalitarian');

MS_d4 = marital_sorting(param,theta_S_d4,Q_d4);
decomp.MS.CMd4 = MS_d4;

avgPWfd4 = sum(sum(MS_d4.*PWf_d4));

[Cf_d4,Cm_d4,Lf_d4,Lm_d4,G_d4,Hf_d4,Hm_d4,H_d4,EQ_d4,...
    Nf_d4,Nm_d4] = per_sol_married(param,wages.f,wages.m,PWf_d4);
            
decomp.tad4.hfmarried = sum(sum(MS_d4.*Hf_d4));
decomp.tad4.nfmarried = sum(sum(MS_d4.*Nf_d4));
decomp.tad4.lfmarried = sum(sum(MS_d4.*Lf_d4));
decomp.tad4.hmmarried = sum(sum(MS_d4.*Hm_d4));
decomp.tad4.nmmarried = sum(sum(MS_d4.*Nm_d4));
decomp.tad4.lmmarried = sum(sum(MS_d4.*Lm_d4));

% Single female time allocation
param.sigma_c = uw_singlefemale(1);
param.sigma_l = uw_singlefemale(2);
param.sigma_g = uw_singlefemale(3);
[~,lfs,~,hfs,~,nfs] = per_sol_singles(param,wages.f);

decomp.tad4.hfsingle = sum(hfs.*param.Pf);
decomp.tad4.nfsingle = sum(nfs.*param.Pf);
decomp.tad4.lfsingle = sum(lfs.*param.Pf);

% Single men time allocation
param.sigma_c = uw_singlemale(1);
param.sigma_l = uw_singlemale(2);
param.sigma_g = uw_singlemale(3);

[~,lms,~,hms,~,nms] = per_sol_singles(param,wages.m);

decomp.tad4.hmsingle = sum(hms.*param.Pm);
decomp.tad4.nmsingle = sum(nms.*param.Pm);
decomp.tad4.lmsingle = sum(lms.*param.Pm);

% Forward decomposition

%% Fifth stage: the sex ratio jumps to 2010 level

param.theta_0 = data.sexratio.sexratio2010;

param.sigma_c = uw_singlefemale(1);
param.sigma_l = uw_singlefemale(2);
param.sigma_g = uw_singlefemale(3);
USfd5 = uflow_singles(param,wages.f);

param.sigma_c = uw_singlemale(1);
param.sigma_l = uw_singlemale(2);
param.sigma_g = uw_singlemale(3);
USmd5 = uflow_singles(param,wages.m);

param.sigma_c = married_sigma_c;
param.sigma_l = married_sigma_l;
param.sigma_g = 1-param.sigma_c-param.sigma_l;

% Compute equilibrium marriage objects
[theta_S_d5,Q_d5,PWf_d5] = sseq(param,wages,USfd5,USmd5,...
    'bargaining','Egalitarian');

MS_d5 = marital_sorting(param,theta_S_d5,Q_d5);
decomp.MS.CMd5 = MS_d5;

avgPWfd5 = sum(sum(MS_d5.*PWf_d5));

[Cf_d5,Cm_d5,Lf_d5,Lm_d5,G_d5,Hf_d5,Hm_d5,H_d5,EQ_d5,...
    Nf_d5,Nm_d5] = per_sol_married(param,wages.f,wages.m,PWf_d5);
            
decomp.tad5.hfmarried = sum(sum(MS_d5.*Hf_d5));
decomp.tad5.nfmarried = sum(sum(MS_d5.*Nf_d5));
decomp.tad5.lfmarried = sum(sum(MS_d5.*Lf_d5));
decomp.tad5.hmmarried = sum(sum(MS_d5.*Hm_d5));
decomp.tad5.nmmarried = sum(sum(MS_d5.*Nm_d5));
decomp.tad5.lmmarried = sum(sum(MS_d5.*Lm_d5));

% Single female time allocation
param.sigma_c = uw_singlefemale(1);
param.sigma_l = uw_singlefemale(2);
param.sigma_g = uw_singlefemale(3);
[~,lfs,~,hfs,~,nfs] = per_sol_singles(param,wages.f);

decomp.tad5.hfsingle = sum(hfs.*param.Pf);
decomp.tad5.nfsingle = sum(nfs.*param.Pf);
decomp.tad5.lfsingle = sum(lfs.*param.Pf);

% Single men time allocation
param.sigma_c = uw_singlemale(1);
param.sigma_l = uw_singlemale(2);
param.sigma_g = uw_singlemale(3);

[~,lms,~,hms,~,nms] = per_sol_singles(param,wages.m);

decomp.tad5.hmsingle = sum(hms.*param.Pm);
decomp.tad5.nmsingle = sum(nms.*param.Pm);
decomp.tad5.lmsingle = sum(lms.*param.Pm);

%% Sixth stage: skill distributions jump to 2010 level

param.Pf = data.skilldist.Pf2010;
param.Pm = data.skilldist.Pm2010;

param.sigma_c = uw_singlefemale(1);
param.sigma_l = uw_singlefemale(2);
param.sigma_g = uw_singlefemale(3);
USfd6 = uflow_singles(param,wages.f);

param.sigma_c = uw_singlemale(1);
param.sigma_l = uw_singlemale(2);
param.sigma_g = uw_singlemale(3);
USmd6 = uflow_singles(param,wages.m);

param.sigma_c = married_sigma_c;
param.sigma_l = married_sigma_l;
param.sigma_g = 1-param.sigma_c-param.sigma_l;

% Compute equilibrium marriage objects
[theta_S_d6,Q_d6,PWf_d6] = sseq(param,wages,USfd6,USmd6,...
    'bargaining','Egalitarian');

MS_d6 = marital_sorting(param,theta_S_d6,Q_d6);
decomp.MS.CMd6 = MS_d6;

avgPWfd6 = sum(sum(MS_d6.*PWf_d6));

[Cf_d6,Cm_d6,Lf_d6,Lm_d6,G_d6,Hf_d6,Hm_d6,H_d6,EQ_d6,...
    Nf_d6,Nm_d6] = per_sol_married(param,wages.f,wages.m,PWf_d6);
            
decomp.tad6.hfmarried = sum(sum(MS_d6.*Hf_d6));
decomp.tad6.nfmarried = sum(sum(MS_d6.*Nf_d6));
decomp.tad6.lfmarried = sum(sum(MS_d6.*Lf_d6));
decomp.tad6.hmmarried = sum(sum(MS_d6.*Hm_d6));
decomp.tad6.nmmarried = sum(sum(MS_d6.*Nm_d6));
decomp.tad6.lmmarried = sum(sum(MS_d6.*Lm_d6));

% Single female time allocation
param.sigma_c = uw_singlefemale(1);
param.sigma_l = uw_singlefemale(2);
param.sigma_g = uw_singlefemale(3);
[~,lfs,~,hfs,~,nfs] = per_sol_singles(param,wages.f);

decomp.tad6.hfsingle = sum(hfs.*param.Pf);
decomp.tad6.nfsingle = sum(nfs.*param.Pf);
decomp.tad6.lfsingle = sum(lfs.*param.Pf);

% Single men time allocation
param.sigma_c = uw_singlemale(1);
param.sigma_l = uw_singlemale(2);
param.sigma_g = uw_singlemale(3);

[~,lms,~,hms,~,nms] = per_sol_singles(param,wages.m);

decomp.tad6.hmsingle = sum(hms.*param.Pm);
decomp.tad6.nmsingle = sum(nms.*param.Pm);
decomp.tad6.lmsingle = sum(lms.*param.Pm);

%% Seventh stage: wages jumo to 2010 level

wages.m = wages.m2010;
wages.f = wages.f2010;

param.sigma_c = uw_singlefemale(1);
param.sigma_l = uw_singlefemale(2);
param.sigma_g = uw_singlefemale(3);
USfd7 = uflow_singles(param,wages.f);

param.sigma_c = uw_singlemale(1);
param.sigma_l = uw_singlemale(2);
param.sigma_g = uw_singlemale(3);
USmd7 = uflow_singles(param,wages.m);

param.sigma_c = married_sigma_c;
param.sigma_l = married_sigma_l;
param.sigma_g = 1-param.sigma_c-param.sigma_l;

% Compute equilibrium marriage objects
[theta_S_d7,Q_d7,PWf_d7] = sseq(param,wages,USfd7,USmd7,...
    'bargaining','Egalitarian');

MS_d7 = marital_sorting(param,theta_S_d7,Q_d7);
decomp.MS.CMd7 = MS_d7;

avgPWfd7 = sum(sum(MS_d7.*PWf_d7));

[Cf_d7,Cm_d7,Lf_d7,Lm_d7,G_d7,Hf_d7,Hm_d7,H_d7,EQ_d7,...
    Nf_d7,Nm_d7] = per_sol_married(param,wages.f,wages.m,PWf_d7);
            
decomp.tad7.hfmarried = sum(sum(MS_d7.*Hf_d7));
decomp.tad7.nfmarried = sum(sum(MS_d7.*Nf_d7));
decomp.tad7.lfmarried = sum(sum(MS_d7.*Lf_d7));
decomp.tad7.hmmarried = sum(sum(MS_d7.*Hm_d7));
decomp.tad7.nmmarried = sum(sum(MS_d7.*Nm_d7));
decomp.tad7.lmmarried = sum(sum(MS_d7.*Lm_d7));

% Single female time allocation
param.sigma_c = uw_singlefemale(1);
param.sigma_l = uw_singlefemale(2);
param.sigma_g = uw_singlefemale(3);
[~,lfs,~,hfs,~,nfs] = per_sol_singles(param,wages.f);

decomp.tad7.hfsingle = sum(hfs.*param.Pf);
decomp.tad7.nfsingle = sum(nfs.*param.Pf);
decomp.tad7.lfsingle = sum(lfs.*param.Pf);

% Single men time allocation
param.sigma_c = uw_singlemale(1);
param.sigma_l = uw_singlemale(2);
param.sigma_g = uw_singlemale(3);

[~,lms,~,hms,~,nms] = per_sol_singles(param,wages.m);

decomp.tad7.hmsingle = sum(hms.*param.Pm);
decomp.tad7.nmsingle = sum(nms.*param.Pm);
decomp.tad7.lmsingle = sum(lms.*param.Pm);

%% Eight stage: pe and Ag jump to 2010 level

param.pe = 0.5;
param.A_g = (1+param.growth_A_g)^20;

param.sigma_c = uw_singlefemale(1);
param.sigma_l = uw_singlefemale(2);
param.sigma_g = uw_singlefemale(3);
USfd8 = uflow_singles(param,wages.f);

param.sigma_c = uw_singlemale(1);
param.sigma_l = uw_singlemale(2);
param.sigma_g = uw_singlemale(3);
USmd8 = uflow_singles(param,wages.m);

param.sigma_c = married_sigma_c;
param.sigma_l = married_sigma_l;
param.sigma_g = 1-param.sigma_c-param.sigma_l;

% Compute equilibrium marriage objects
[theta_S_d8,Q_d8,PWf_d8] = sseq(param,wages,USfd8,USmd8,...
    'bargaining','Egalitarian');

MS_d8 = marital_sorting(param,theta_S_d8,Q_d8);
decomp.MS.CMd8 = MS_d8;

avgPWfd8 = sum(sum(MS_d8.*PWf_d8));

[Cf_d8,Cm_d8,Lf_d8,Lm_d8,G_d8,Hf_d8,Hm_d8,H_d8,EQ_d8,...
    Nf_d8,Nm_d8] = per_sol_married(param,wages.f,wages.m,PWf_d8);
            
decomp.tad8.hfmarried = sum(sum(MS_d8.*Hf_d8));
decomp.tad8.nfmarried = sum(sum(MS_d8.*Nf_d8));
decomp.tad8.lfmarried = sum(sum(MS_d8.*Lf_d8));
decomp.tad8.hmmarried = sum(sum(MS_d8.*Hm_d8));
decomp.tad8.nmmarried = sum(sum(MS_d8.*Nm_d8));
decomp.tad8.lmmarried = sum(sum(MS_d8.*Lm_d8));

% Single female time allocation
param.sigma_c = uw_singlefemale(1);
param.sigma_l = uw_singlefemale(2);
param.sigma_g = uw_singlefemale(3);
[~,lfs,~,hfs,~,nfs] = per_sol_singles(param,wages.f);

decomp.tad8.hfsingle = sum(hfs.*param.Pf);
decomp.tad8.nfsingle = sum(nfs.*param.Pf);
decomp.tad8.lfsingle = sum(lfs.*param.Pf);

% Single men time allocation
param.sigma_c = uw_singlemale(1);
param.sigma_l = uw_singlemale(2);
param.sigma_g = uw_singlemale(3);

[~,lms,~,hms,~,nms] = per_sol_singles(param,wages.m);

decomp.tad8.hmsingle = sum(hms.*param.Pm);
decomp.tad8.nmsingle = sum(nms.*param.Pm);
decomp.tad8.lmsingle = sum(lms.*param.Pm);

%% Exporting the decomposition results

cd(main_dir);

decomp_res.Statistic = {'Married women housework';...
    'Married women paid work';'Married women leisure';...
    'Married men housework';'Married men paid work';...
    'Married men leisure';'Single women housework';...
    'Single women paid work';'Single women leisure';...
    'Single men housework';'Single men paid work';...
    'Single men leisure';'Married women consumption';...
    'Married men consumption';'Average wife Pareto weight';...
    'Assortative mating measure'};

decomp_res.Model2010 = [struct2array(model.ta2010)'*118;...
    sum(Cf_2010.*MS_2010,'all');sum(Cm_2010.*MS_2010,'all');...
    avgPWf2010;msmeasures(model.MS.CM2010)];

decomp_res.D1 = [struct2array(decomp.tad1)'*118;...
    sum(Cf_d1.*MS_d1,'all');sum(Cm_d1.*MS_d1,'all');...
    avgPWfd1;msmeasures(decomp.MS.CMd1)];

decomp_res.D2 = [struct2array(decomp.tad2)'*118;...
    sum(Cf_d2.*MS_d2,'all');sum(Cm_d2.*MS_d2,'all');...
    avgPWfd2;msmeasures(decomp.MS.CMd2)];

decomp_res.D3 = [struct2array(decomp.tad3)'*118;...
    sum(Cf_d3.*MS_d3,'all');sum(Cm_d3.*MS_d3,'all');...
    avgPWfd3;msmeasures(decomp.MS.CMd3)];

decomp_res.D4 = [struct2array(decomp.tad4)'*118;...
    sum(Cf_d4.*MS_d4,'all');sum(Cm_d4.*MS_d4,'all');...
    avgPWfd4;msmeasures(decomp.MS.CMd4)];

decomp_res.D5 = [struct2array(decomp.tad5)'*118;...
    sum(Cf_d5.*MS_d5,'all');sum(Cm_d5.*MS_d5,'all');...
    avgPWfd5;msmeasures(decomp.MS.CMd5)];

decomp_res.D6 = [struct2array(decomp.tad6)'*118;...
    sum(Cf_d6.*MS_d6,'all');sum(Cm_d6.*MS_d6,'all');...
    avgPWfd6;msmeasures(decomp.MS.CMd6)];

decomp_res.D7 = [struct2array(decomp.tad7)'*118;...
    sum(Cf_d7.*MS_d7,'all');sum(Cm_d7.*MS_d7,'all');...
    avgPWfd7;msmeasures(decomp.MS.CMd7)];

decomp_res.D8 = [struct2array(decomp.tad8)'*118;...
    sum(Cf_d8.*MS_d8,'all');sum(Cm_d8.*MS_d8,'all');...
    avgPWfd8;msmeasures(decomp.MS.CMd8)];

decomposition_results = struct2table(decomp_res);

cd(results_dir);

Sheet = "decomp_res, sigma=" + param.sigma;

writetable(decomposition_results,filename_results,...
    'Sheet',Sheet);

%% Figures for the decomposition

cd(results_dir);
Sheet = "decomp_res_perc, sigma=" + param.sigma;

decomp_el = {'Sex ratio';'Skill distribution';'Wages';'Home production'};

bar_color = [0.75 0.75 0.75];
line_color = 'k';

decomp_back = readmatrix(filename_results,'Sheet',Sheet,'Range','B2:F17');
decomp_for = readmatrix(filename_results,'Sheet',Sheet,'Range','B19:F34');

%% Backward decomposition time allocation females

housework_married_female = decomp_back(1,:)*100;
paid_work_married_female = decomp_back(2,:)*100;
leisure_married_female = decomp_back(3,:)*100;

decomp_back_diff = decomp_back(:,2:end)-decomp_back(:,1:end-1);

housework_married_female_diff = decomp_back_diff(1,:)*100;
paid_work_married_female_diff = decomp_back_diff(2,:)*100;
leisure_married_female_diff = decomp_back_diff(3,:)*100;

txt_sexratio_housework = [num2str(housework_married_female_diff(1),'%.0f') '%'];
txt_skills_housework = [num2str(housework_married_female_diff(2),'%.0f') '%'];
txt_wages_housework = [num2str(housework_married_female_diff(3),'%.0f') '%'];
txt_home_housework = [num2str(housework_married_female_diff(4),'%.0f') '%'];

txt_sexratio_paid_work = [num2str(paid_work_married_female_diff(1),'%.0f') '%'];
txt_skills_paid_work = [num2str(paid_work_married_female_diff(2),'%.0f') '%'];
txt_wages_paid_work = [num2str(paid_work_married_female_diff(3),'%.0f') '%'];
txt_home_paid_work = [num2str(paid_work_married_female_diff(4),'%.0f') '%'];

txt_sexratio_leisure = [num2str(leisure_married_female_diff(1),'%.0f') '%'];
txt_skills_leisure = [num2str(leisure_married_female_diff(2),'%.0f') '%'];
txt_wages_leisure = [num2str(leisure_married_female_diff(3),'%.0f') '%'];
txt_home_leisure = [num2str(leisure_married_female_diff(4),'%.0f') '%'];

figure('Color','White');
tiledlayout(3,1);

% Top plot: married women housework
ax1 = nexttile;
hold on 
bar(2,housework_married_female_diff(1),'FaceColor',bar_color)
bar(3,housework_married_female_diff(2),'FaceColor',bar_color)
bar(4,housework_married_female_diff(3),'FaceColor',bar_color)
bar(5,housework_married_female_diff(4),'FaceColor',bar_color)
plot(housework_married_female,'-o','Color',line_color,'MarkerFaceColor',line_color)
plot([2,2],[-100,200],'--','Color','red')
plot([3,3],[-100,200],'--','Color','red')
plot([4,4],[-100,200],'--','Color','red')
plot([5,5],[-100,200],'--','Color','red')
yline(100,'red')
text(2,housework_married_female_diff(1),txt_sexratio_housework,...
    'HorizontalAlignment','center','VerticalAlignment', 'bottom')
text(3,housework_married_female_diff(2),txt_skills_housework,...
    'HorizontalAlignment','center','VerticalAlignment', 'top')
text(4,housework_married_female_diff(3),txt_wages_housework,...
    'HorizontalAlignment','center','VerticalAlignment', 'bottom')
text(5,housework_married_female_diff(4),txt_home_housework,...
    'HorizontalAlignment','center','VerticalAlignment', 'bottom')
set(gca,'Color','white','xtick',[2,3,4,5],'xticklabel',decomp_el,'TickLength',[0 0])
ytickformat('percentage')
xlim([1 6])
hold off
title(ax1,'Housework')

% Middle plot: married women paid work
ax2 = nexttile;
hold on 
bar(2,paid_work_married_female_diff(1),'FaceColor',bar_color)
bar(3,paid_work_married_female_diff(2),'FaceColor',bar_color)
bar(4,paid_work_married_female_diff(3),'FaceColor',bar_color)
bar(5,paid_work_married_female_diff(4),'FaceColor',bar_color)
plot(paid_work_married_female,'-o','Color',line_color,'MarkerFaceColor',line_color)
plot([2,2],[-100,200],'--','Color','red')
plot([3,3],[-100,200],'--','Color','red')
plot([4,4],[-100,200],'--','Color','red')
plot([5,5],[-100,200],'--','Color','red')
yline(100,'red')
text(2,paid_work_married_female_diff(1),txt_sexratio_paid_work,...
    'HorizontalAlignment','center','VerticalAlignment', 'bottom')
text(3,paid_work_married_female_diff(2),txt_skills_paid_work,...
    'HorizontalAlignment','center','VerticalAlignment', 'top')
text(4,paid_work_married_female_diff(3),txt_wages_paid_work,...
    'HorizontalAlignment','center','VerticalAlignment', 'top')
text(5,paid_work_married_female_diff(4),txt_home_paid_work,...
    'HorizontalAlignment','center','VerticalAlignment', 'bottom')
set(gca,'Color','white','xtick',[2,3,4,5],'xticklabel',decomp_el,'TickLength',[0 0])
ytickformat('percentage')
xlim([1 6])
hold off
title(ax2,'Paid work')

% Bottom plot: married women leisure
ax3 = nexttile;
hold on
bar(2,leisure_married_female_diff(1),'FaceColor',bar_color)
bar(3,leisure_married_female_diff(2),'FaceColor',bar_color)
bar(4,leisure_married_female_diff(3),'FaceColor',bar_color)
bar(5,leisure_married_female_diff(4),'FaceColor',bar_color)
plot(leisure_married_female,'-o','Color',line_color,'MarkerFaceColor',line_color)
plot([2,2],[-200,100],'--','Color','red')
plot([3,3],[-200,100],'--','Color','red')
plot([4,4],[-200,100],'--','Color','red')
plot([5,5],[-200,100],'--','Color','red')
yline(-100,'red')
text(2,leisure_married_female_diff(1),txt_sexratio_leisure,...
    'HorizontalAlignment','center','VerticalAlignment', 'top')
text(3,leisure_married_female_diff(2),txt_skills_leisure,...
    'HorizontalAlignment','center','VerticalAlignment', 'bottom')
text(4,leisure_married_female_diff(3),txt_wages_leisure,...
    'HorizontalAlignment','center','VerticalAlignment', 'top')
text(5,leisure_married_female_diff(4),txt_home_leisure,...
    'HorizontalAlignment','center','VerticalAlignment', 'top')
set(gca,'Color','white','xtick',[2,3,4,5],'xticklabel',decomp_el,'TickLength',[0 0])
ytickformat('percentage')
xlim([1 6])
hold off
title(ax3,'Leisure')

saveas(gcf,'back_decomp_ta_f','epsc')

%% Forward decomposition time allocation females

housework_married_female = decomp_for(1,:)*100;
paid_work_married_female = decomp_for(2,:)*100;
leisure_married_female = decomp_for(3,:)*100;

decomp_for_diff = decomp_for(:,2:end)-decomp_for(:,1:end-1);

housework_married_female_diff = decomp_for_diff(1,:)*100;
paid_work_married_female_diff = decomp_for_diff(2,:)*100;
leisure_married_female_diff = decomp_for_diff(3,:)*100;

txt_sexratio_housework = [num2str(housework_married_female_diff(1),'%.0f') '%'];
txt_skills_housework = [num2str(housework_married_female_diff(2),'%.0f') '%'];
txt_wages_housework = [num2str(housework_married_female_diff(3),'%.0f') '%'];
txt_home_housework = [num2str(housework_married_female_diff(4),'%.0f') '%'];

txt_sexratio_paid_work = [num2str(paid_work_married_female_diff(1),'%.0f') '%'];
txt_skills_paid_work = [num2str(paid_work_married_female_diff(2),'%.0f') '%'];
txt_wages_paid_work = [num2str(paid_work_married_female_diff(3),'%.0f') '%'];
txt_home_paid_work = [num2str(paid_work_married_female_diff(4),'%.0f') '%'];

txt_sexratio_leisure = [num2str(leisure_married_female_diff(1),'%.0f') '%'];
txt_skills_leisure = [num2str(leisure_married_female_diff(2),'%.0f') '%'];
txt_wages_leisure = [num2str(leisure_married_female_diff(3),'%.0f') '%'];
txt_home_leisure = [num2str(leisure_married_female_diff(4),'%.0f') '%'];

figure('Color','White');
tiledlayout(3,1);

% Top plot: married women housework
ax1 = nexttile;
hold on 
bar(2,housework_married_female_diff(1),'FaceColor',bar_color)
bar(3,housework_married_female_diff(2),'FaceColor',bar_color)
bar(4,housework_married_female_diff(3),'FaceColor',bar_color)
bar(5,housework_married_female_diff(4),'FaceColor',bar_color)
plot(housework_married_female,'-o','Color',line_color,'MarkerFaceColor',line_color)
plot([2,2],[-200,200],'--','Color','red')
plot([3,3],[-200,200],'--','Color','red')
plot([4,4],[-200,200],'--','Color','red')
plot([5,5],[-200,200],'--','Color','red')
yline(-100,'red')
text(2,housework_married_female_diff(1),txt_sexratio_housework,...
    'HorizontalAlignment','center','VerticalAlignment', 'bottom')
text(3,housework_married_female_diff(2),txt_skills_housework,...
    'HorizontalAlignment','center','VerticalAlignment', 'top')
text(4,housework_married_female_diff(3),txt_wages_housework,...
    'HorizontalAlignment','center','VerticalAlignment', 'bottom')
text(5,housework_married_female_diff(4),txt_home_housework,...
    'HorizontalAlignment','center','VerticalAlignment', 'bottom')
set(gca,'Color','white','xtick',[2,3,4,5],'xticklabel',decomp_el,'TickLength',[0 0])
ytickformat('percentage')
xlim([1 6])
hold off
title(ax1,'Housework')

% Middle plot: married women paid work
ax2 = nexttile;
hold on 
bar(2,paid_work_married_female_diff(1),'FaceColor',bar_color)
bar(3,paid_work_married_female_diff(2),'FaceColor',bar_color)
bar(4,paid_work_married_female_diff(3),'FaceColor',bar_color)
bar(5,paid_work_married_female_diff(4),'FaceColor',bar_color)
plot(paid_work_married_female,'-o','Color',line_color,'MarkerFaceColor',line_color)
plot([2,2],[-200,100],'--','Color','red')
plot([3,3],[-200,100],'--','Color','red')
plot([4,4],[-200,100],'--','Color','red')
plot([5,5],[-200,100],'--','Color','red')
yline(-100,'red')
text(2,paid_work_married_female_diff(1),txt_sexratio_paid_work,...
    'HorizontalAlignment','center','VerticalAlignment', 'top')
text(3,paid_work_married_female_diff(2),txt_skills_paid_work,...
    'HorizontalAlignment','center','VerticalAlignment', 'bottom')
text(4,paid_work_married_female_diff(3),txt_wages_paid_work,...
    'HorizontalAlignment','center','VerticalAlignment', 'bottom')
text(5,paid_work_married_female_diff(4),txt_home_paid_work,...
    'HorizontalAlignment','center','VerticalAlignment', 'top')
set(gca,'Color','white','xtick',[2,3,4,5],'xticklabel',decomp_el,'TickLength',[0 0])
ytickformat('percentage')
xlim([1 6])
hold off
title(ax2,'Paid work')

% Bottom plot: married women leisure
ax3 = nexttile;
hold on
bar(2,leisure_married_female_diff(1),'FaceColor',bar_color)
bar(3,leisure_married_female_diff(2),'FaceColor',bar_color)
bar(4,leisure_married_female_diff(3),'FaceColor',bar_color)
bar(5,leisure_married_female_diff(4),'FaceColor',bar_color)
plot(leisure_married_female,'-o','Color',line_color,'MarkerFaceColor',line_color)
plot([2,2],[-100,200],'--','Color','red')
plot([3,3],[-100,200],'--','Color','red')
plot([4,4],[-100,200],'--','Color','red')
plot([5,5],[-100,200],'--','Color','red')
yline(100,'red')
text(2,leisure_married_female_diff(1),txt_sexratio_leisure,...
    'HorizontalAlignment','center','VerticalAlignment', 'bottom')
text(3,leisure_married_female_diff(2),txt_skills_leisure,...
    'HorizontalAlignment','center','VerticalAlignment', 'top')
text(4,leisure_married_female_diff(3),txt_wages_leisure,...
    'HorizontalAlignment','center','VerticalAlignment', 'top')
text(5,leisure_married_female_diff(4),txt_home_leisure,...
    'HorizontalAlignment','center','VerticalAlignment', 'bottom')
set(gca,'Color','white','xtick',[2,3,4,5],'xticklabel',decomp_el,'TickLength',[0 0])
ytickformat('percentage')
xlim([1 6])
hold off
title(ax3,'Leisure')

saveas(gcf,'for_decomp_ta_f','epsc')

%% Backward decomposition time allocation males

housework_married_male = decomp_back(4,:)*100;
paid_work_married_male = decomp_back(5,:)*100;
leisure_married_male = decomp_back(6,:)*100;

housework_married_male_diff = decomp_back_diff(4,:)*100;
paid_work_married_male_diff = decomp_back_diff(5,:)*100;
leisure_married_male_diff = decomp_back_diff(6,:)*100;

txt_sexratio_housework = [num2str(housework_married_male_diff(1),'%.0f') '%'];
txt_skills_housework = [num2str(housework_married_male_diff(2),'%.0f') '%'];
txt_wages_housework = [num2str(housework_married_male_diff(3),'%.0f') '%'];
txt_home_housework = [num2str(housework_married_male_diff(4),'%.0f') '%'];

txt_sexratio_paid_work = [num2str(paid_work_married_male_diff(1),'%.0f') '%'];
txt_skills_paid_work = [num2str(paid_work_married_male_diff(2),'%.0f') '%'];
txt_wages_paid_work = [num2str(paid_work_married_male_diff(3),'%.0f') '%'];
txt_home_paid_work = [num2str(paid_work_married_male_diff(4),'%.0f') '%'];

txt_sexratio_leisure = [num2str(leisure_married_male_diff(1),'%.0f') '%'];
txt_skills_leisure = [num2str(leisure_married_male_diff(2),'%.0f') '%'];
txt_wages_leisure = [num2str(leisure_married_male_diff(3),'%.0f') '%'];
txt_home_leisure = [num2str(leisure_married_male_diff(4),'%.0f') '%'];

figure('Color','White');
tiledlayout(3,1);

% Top plot: married men housework
ax1 = nexttile;
hold on 
bar(2,housework_married_male_diff(1),'FaceColor',bar_color)
bar(3,housework_married_male_diff(2),'FaceColor',bar_color)
bar(4,housework_married_male_diff(3),'FaceColor',bar_color)
bar(5,housework_married_male_diff(4),'FaceColor',bar_color)
plot(housework_married_male,'-o','Color',line_color,'MarkerFaceColor',line_color)
plot([2,2],[-100,200],'--','Color','red')
plot([3,3],[-100,200],'--','Color','red')
plot([4,4],[-100,200],'--','Color','red')
plot([5,5],[-100,200],'--','Color','red')
yline(100,'red')
text(2,housework_married_male_diff(1),txt_sexratio_housework,...
    'HorizontalAlignment','center','VerticalAlignment', 'top')
text(3,housework_married_male_diff(2),txt_skills_housework,...
    'HorizontalAlignment','center','VerticalAlignment', 'top')
text(4,housework_married_male_diff(3),txt_wages_housework,...
    'HorizontalAlignment','center','VerticalAlignment', 'bottom')
text(5,housework_married_male_diff(4),txt_home_housework,...
    'HorizontalAlignment','center','VerticalAlignment', 'bottom')
set(gca,'Color','white','xtick',[2,3,4,5],'xticklabel',decomp_el,'TickLength',[0 0])
ytickformat('percentage')
xlim([1 6])
hold off
title(ax1,'Housework')

% Middle plot: married men paid work
ax2 = nexttile;
hold on 
bar(2,paid_work_married_male_diff(1),'FaceColor',bar_color)
bar(3,paid_work_married_male_diff(2),'FaceColor',bar_color)
bar(4,paid_work_married_male_diff(3),'FaceColor',bar_color)
bar(5,paid_work_married_male_diff(4),'FaceColor',bar_color)
plot(paid_work_married_male,'-o','Color',line_color,'MarkerFaceColor',line_color)
plot([2,2],[-200,200],'--','Color','red')
plot([3,3],[-200,200],'--','Color','red')
plot([4,4],[-200,200],'--','Color','red')
plot([5,5],[-200,200],'--','Color','red')
yline(100,'red')
text(2,paid_work_married_male_diff(1),txt_sexratio_paid_work,...
    'HorizontalAlignment','center','VerticalAlignment', 'top')
text(3,paid_work_married_male_diff(2),txt_skills_paid_work,...
    'HorizontalAlignment','center','VerticalAlignment', 'bottom')
text(4,paid_work_married_male_diff(3),txt_wages_paid_work,...
    'HorizontalAlignment','center','VerticalAlignment', 'top')
text(5,paid_work_married_male_diff(4),txt_home_paid_work,...
    'HorizontalAlignment','center','VerticalAlignment', 'top')
set(gca,'Color','white','xtick',[2,3,4,5],'xticklabel',decomp_el,'TickLength',[0 0])
ytickformat('percentage')
xlim([1 6])
hold off
title(ax2,'Paid work')

% Bottom plot: married men leisure
ax3 = nexttile;
hold on
bar(2,leisure_married_male_diff(1),'FaceColor',bar_color)
bar(3,leisure_married_male_diff(2),'FaceColor',bar_color)
bar(4,leisure_married_male_diff(3),'FaceColor',bar_color)
bar(5,leisure_married_male_diff(4),'FaceColor',bar_color)
plot(leisure_married_male,'-o','Color',line_color,'MarkerFaceColor',line_color)
plot([2,2],[-200,200],'--','Color','red')
plot([3,3],[-200,200],'--','Color','red')
plot([4,4],[-200,200],'--','Color','red')
plot([5,5],[-200,200],'--','Color','red')
yline(-100,'red')
text(2,leisure_married_male_diff(1),txt_sexratio_leisure,...
    'HorizontalAlignment','center','VerticalAlignment', 'bottom')
text(3,leisure_married_male_diff(2),txt_skills_leisure,...
    'HorizontalAlignment','center','VerticalAlignment', 'top')
text(4,leisure_married_male_diff(3),txt_wages_leisure,...
    'HorizontalAlignment','center','VerticalAlignment', 'top')
text(5,leisure_married_male_diff(4),txt_home_leisure,...
    'HorizontalAlignment','center','VerticalAlignment', 'bottom')
set(gca,'Color','white','xtick',[2,3,4,5],'xticklabel',decomp_el,'TickLength',[0 0])
ytickformat('percentage')
xlim([1 6])
hold off
title(ax3,'Leisure')

saveas(gcf,'back_decomp_ta_m','epsc')

%% Forward decomposition time allocation males

housework_married_male = decomp_for(4,:)*100;
paid_work_married_male = decomp_for(5,:)*100;
leisure_married_male = decomp_for(6,:)*100;

housework_married_male_diff = decomp_for_diff(4,:)*100;
paid_work_married_male_diff = decomp_for_diff(5,:)*100;
leisure_married_male_diff = decomp_for_diff(6,:)*100;

txt_sexratio_housework = [num2str(housework_married_male_diff(1),'%.0f') '%'];
txt_skills_housework = [num2str(housework_married_male_diff(2),'%.0f') '%'];
txt_wages_housework = [num2str(housework_married_male_diff(3),'%.0f') '%'];
txt_home_housework = [num2str(housework_married_male_diff(4),'%.0f') '%'];

txt_sexratio_paid_work = [num2str(paid_work_married_male_diff(1),'%.0f') '%'];
txt_skills_paid_work = [num2str(paid_work_married_male_diff(2),'%.0f') '%'];
txt_wages_paid_work = [num2str(paid_work_married_male_diff(3),'%.0f') '%'];
txt_home_paid_work = [num2str(paid_work_married_male_diff(4),'%.0f') '%'];

txt_sexratio_leisure = [num2str(leisure_married_male_diff(1),'%.0f') '%'];
txt_skills_leisure = [num2str(leisure_married_male_diff(2),'%.0f') '%'];
txt_wages_leisure = [num2str(leisure_married_male_diff(3),'%.0f') '%'];
txt_home_leisure = [num2str(leisure_married_male_diff(4),'%.0f') '%'];

figure('Color','White');
tiledlayout(3,1);

% Top plot: married men housework
ax1 = nexttile;
hold on 
bar(2,housework_married_male_diff(1),'FaceColor',bar_color)
bar(3,housework_married_male_diff(2),'FaceColor',bar_color)
bar(4,housework_married_male_diff(3),'FaceColor',bar_color)
bar(5,housework_married_male_diff(4),'FaceColor',bar_color)
plot(housework_married_male,'-o','Color',line_color,'MarkerFaceColor',line_color)
plot([2,2],[-200,100],'--','Color','red')
plot([3,3],[-200,100],'--','Color','red')
plot([4,4],[-200,100],'--','Color','red')
plot([5,5],[-200,100],'--','Color','red')
yline(-100,'red')
text(2,housework_married_male_diff(1),txt_sexratio_housework,...
    'HorizontalAlignment','center','VerticalAlignment', 'bottom')
text(3,housework_married_male_diff(2),txt_skills_housework,...
    'HorizontalAlignment','center','VerticalAlignment', 'bottom')
text(4,housework_married_male_diff(3),txt_wages_housework,...
    'HorizontalAlignment','center','VerticalAlignment', 'top')
text(5,housework_married_male_diff(4),txt_home_housework,...
    'HorizontalAlignment','center','VerticalAlignment', 'top')
set(gca,'Color','white','xtick',[2,3,4,5],'xticklabel',decomp_el,'TickLength',[0 0])
ytickformat('percentage')
xlim([1 6])
hold off
title(ax1,'Housework')

% Middle plot: married men paid work
ax2 = nexttile;
hold on 
bar(2,paid_work_married_male_diff(1),'FaceColor',bar_color)
bar(3,paid_work_married_male_diff(2),'FaceColor',bar_color)
bar(4,paid_work_married_male_diff(3),'FaceColor',bar_color)
bar(5,paid_work_married_male_diff(4),'FaceColor',bar_color)
plot(paid_work_married_male,'-o','Color',line_color,'MarkerFaceColor',line_color)
plot([2,2],[-300,300],'--','Color','red')
plot([3,3],[-300,300],'--','Color','red')
plot([4,4],[-300,300],'--','Color','red')
plot([5,5],[-300,300],'--','Color','red')
yline(-100,'red')
text(2,paid_work_married_male_diff(1),txt_sexratio_paid_work,...
    'HorizontalAlignment','center','VerticalAlignment', 'bottom')
text(3,paid_work_married_male_diff(2),txt_skills_paid_work,...
    'HorizontalAlignment','center','VerticalAlignment', 'top')
text(4,paid_work_married_male_diff(3),txt_wages_paid_work,...
    'HorizontalAlignment','center','VerticalAlignment', 'top')
text(5,paid_work_married_male_diff(4),txt_home_paid_work,...
    'HorizontalAlignment','center','VerticalAlignment', 'bottom')
set(gca,'Color','white','xtick',[2,3,4,5],'xticklabel',decomp_el,'TickLength',[0 0])
ytickformat('percentage')
xlim([1 6])
hold off
title(ax2,'Paid work')

% Bottom plot: married men leisure
ax3 = nexttile;
hold on
bar(2,leisure_married_male_diff(1),'FaceColor',bar_color)
bar(3,leisure_married_male_diff(2),'FaceColor',bar_color)
bar(4,leisure_married_male_diff(3),'FaceColor',bar_color)
bar(5,leisure_married_male_diff(4),'FaceColor',bar_color)
plot(leisure_married_male,'-o','Color',line_color,'MarkerFaceColor',line_color)
plot([2,2],[-200,200],'--','Color','red')
plot([3,3],[-200,200],'--','Color','red')
plot([4,4],[-200,200],'--','Color','red')
plot([5,5],[-200,200],'--','Color','red')
yline(100,'red')
text(2,leisure_married_male_diff(1),txt_sexratio_leisure,...
    'HorizontalAlignment','center','VerticalAlignment', 'top')
text(3,leisure_married_male_diff(2),txt_skills_leisure,...
    'HorizontalAlignment','center','VerticalAlignment', 'bottom')
text(4,leisure_married_male_diff(3),txt_wages_leisure,...
    'HorizontalAlignment','center','VerticalAlignment', 'bottom')
text(5,leisure_married_male_diff(4),txt_home_leisure,...
    'HorizontalAlignment','center','VerticalAlignment', 'top')
set(gca,'Color','white','xtick',[2,3,4,5],'xticklabel',decomp_el,'TickLength',[0 0])
ytickformat('percentage')
xlim([1 6])
hold off
title(ax3,'Leisure')

saveas(gcf,'for_decomp_ta_m','epsc')

%% Backward decomposition assortative mating

assortative_mating = decomp_back(16,:)*100;

assortative_mating_diff = decomp_back_diff(16,:)*100;

txt_sexratio_assortative_mating = [num2str(assortative_mating_diff(1),'%.0f') '%'];
txt_skills_assortative_mating = [num2str(assortative_mating_diff(2),'%.0f') '%'];
txt_wages_assortative_mating = [num2str(assortative_mating_diff(3),'%.0f') '%'];
txt_home_assortative_mating = [num2str(assortative_mating_diff(4),'%.0f') '%'];

figure('Color','White');

hold on 
bar(2,assortative_mating_diff(1),'FaceColor',bar_color)
bar(3,assortative_mating_diff(2),'FaceColor',bar_color)
bar(4,assortative_mating_diff(3),'FaceColor',bar_color)
bar(5,assortative_mating_diff(4),'FaceColor',bar_color)
plot(assortative_mating,'-o','Color',line_color,'MarkerFaceColor',line_color)
plot([2,2],[-200,100],'--','Color','red')
plot([3,3],[-200,100],'--','Color','red')
plot([4,4],[-200,100],'--','Color','red')
plot([5,5],[-200,100],'--','Color','red')
yline(-100,'red')
text(2,assortative_mating_diff(1),txt_sexratio_assortative_mating,...
    'HorizontalAlignment','center','VerticalAlignment', 'top')
text(3,assortative_mating_diff(2),txt_skills_assortative_mating,...
    'HorizontalAlignment','center','VerticalAlignment', 'top')
text(4,assortative_mating_diff(3),txt_wages_assortative_mating,...
    'HorizontalAlignment','center','VerticalAlignment', 'bottom')
text(5,assortative_mating_diff(4),txt_home_assortative_mating,...
    'HorizontalAlignment','center','VerticalAlignment', 'bottom')
set(gca,'Color','white','xtick',[2,3,4,5],'xticklabel',decomp_el,'TickLength',[0 0])
ytickformat('percentage')
xlim([1 6])
hold off

saveas(gcf,'back_decomp_ms','epsc')

%% Forward decomposition assortative mating

assortative_mating = decomp_for(16,:)*100;

assortative_mating_diff = decomp_for_diff(16,:)*100;

txt_sexratio_assortative_mating = [num2str(assortative_mating_diff(1),'%.0f') '%'];
txt_skills_assortative_mating = [num2str(assortative_mating_diff(2),'%.0f') '%'];
txt_wages_assortative_mating = [num2str(assortative_mating_diff(3),'%.0f') '%'];
txt_home_assortative_mating = [num2str(assortative_mating_diff(4),'%.0f') '%'];

figure('Color','White');

hold on 
bar(2,assortative_mating_diff(1),'FaceColor',bar_color)
bar(3,assortative_mating_diff(2),'FaceColor',bar_color)
bar(4,assortative_mating_diff(3),'FaceColor',bar_color)
bar(5,assortative_mating_diff(4),'FaceColor',bar_color)
plot(assortative_mating,'-o','Color',line_color,'MarkerFaceColor',line_color)
plot([2,2],[-100,200],'--','Color','red')
plot([3,3],[-100,200],'--','Color','red')
plot([4,4],[-100,200],'--','Color','red')
plot([5,5],[-100,200],'--','Color','red')
yline(100,'red')
text(2,assortative_mating_diff(1),txt_sexratio_assortative_mating,...
    'HorizontalAlignment','center','VerticalAlignment', 'bottom')
text(3,assortative_mating_diff(2),txt_skills_assortative_mating,...
    'HorizontalAlignment','center','VerticalAlignment', 'bottom')
text(4,assortative_mating_diff(3),txt_wages_assortative_mating,...
    'HorizontalAlignment','center','VerticalAlignment', 'top')
text(5,assortative_mating_diff(4),txt_home_assortative_mating,...
    'HorizontalAlignment','center','VerticalAlignment', 'top')
set(gca,'Color','white','xtick',[2,3,4,5],'xticklabel',decomp_el,'TickLength',[0 0])
ytickformat('percentage')
xlim([1 6])
hold off

saveas(gcf,'for_decomp_ms','epsc')

%%  Decomposing the effect of the sex ratio on married people time allocation

cd(main_dir);

param.sigma_c = married_sigma_c;
param.sigma_l = married_sigma_l;
param.sigma_g = 1-param.sigma_c-param.sigma_l;

% Backward

param.pe = 0.5;
param.A_g = (1+param.growth_A_g)^20;

param.theta_0 = data.sexratio.sexratio2010;

param.Pf = data.skilldist.Pf2010;
param.Pm = data.skilldist.Pm2010;

wages.f = wages.f2010;
wages.m = wages.m2010;

[Cf_d11,Cm_d11,Lf_d11,Lm_d11,G_d11,Hf_d11,Hm_d11,H_d11,EQ_d11,...
    Nf_d11,Nm_d11] = per_sol_married(param,wages.f,wages.m,PWf_d1);
            
decomp.tad11.hfmarried = sum(sum(MS_2010.*Hf_d11));
decomp.tad11.nfmarried = sum(sum(MS_2010.*Nf_d11));
decomp.tad11.lfmarried = sum(sum(MS_2010.*Lf_d11));
decomp.tad11.hmmarried = sum(sum(MS_2010.*Hm_d11));
decomp.tad11.nmmarried = sum(sum(MS_2010.*Nm_d11));
decomp.tad11.lmmarried = sum(sum(MS_2010.*Lm_d11));

[Cf_d12,Cm_d12,Lf_d12,Lm_d12,G_d12,Hf_d12,Hm_d12,H_d12,EQ_d12,...
    Nf_d12,Nm_d12] = per_sol_married(param,wages.f,wages.m,PWf_2010);
            
decomp.tad12.hfmarried = sum(sum(MS_d1.*Hf_d12));
decomp.tad12.nfmarried = sum(sum(MS_d1.*Nf_d12));
decomp.tad12.lfmarried = sum(sum(MS_d1.*Lf_d12));
decomp.tad12.hmmarried = sum(sum(MS_d1.*Hm_d12));
decomp.tad12.nmmarried = sum(sum(MS_d1.*Nm_d12));
decomp.tad12.lmmarried = sum(sum(MS_d1.*Lm_d12));

% Forward

param.pe = 1.5;
param.A_g = 1;

param.theta_0 = data.sexratio.sexratio1990;

param.Pf = data.skilldist.Pf1990;
param.Pm = data.skilldist.Pm1990;

wages.m = wages.m1990;
wages.f = wages.f1990;

[Cf_d51,Cm_d51,Lf_d51,Lm_d51,G_d51,Hf_d51,Hm_d51,H_d51,EQ_d51,...
    Nf_d51,Nm_d51] = per_sol_married(param,wages.f,wages.m,PWf_d5);
            
decomp.tad51.hfmarried = sum(sum(MS_1990.*Hf_d51));
decomp.tad51.nfmarried = sum(sum(MS_1990.*Nf_d51));
decomp.tad51.lfmarried = sum(sum(MS_1990.*Lf_d51));
decomp.tad51.hmmarried = sum(sum(MS_1990.*Hm_d51));
decomp.tad51.nmmarried = sum(sum(MS_1990.*Nm_d51));
decomp.tad51.lmmarried = sum(sum(MS_1990.*Lm_d51));

[Cf_d52,Cm_d52,Lf_d52,Lm_d52,G_d52,Hf_d52,Hm_d52,H_d52,EQ_d52,...
    Nf_d52,Nm_d52] = per_sol_married(param,wages.f,wages.m,PWf_1990);
            
decomp.tad52.hfmarried = sum(sum(MS_d5.*Hf_d52));
decomp.tad52.nfmarried = sum(sum(MS_d5.*Nf_d52));
decomp.tad52.lfmarried = sum(sum(MS_d5.*Lf_d52));
decomp.tad52.hmmarried = sum(sum(MS_d5.*Hm_d52));
decomp.tad52.nmmarried = sum(sum(MS_d5.*Nm_d52));
decomp.tad52.lmmarried = sum(sum(MS_d5.*Lm_d52));

% Export the results

decomp_sexratio_for.Statistic = {'Married women housework';...
    'Married women paid work';'Married women leisure';...
    'Married men housework';'Married men paid work';...
    'Married men leisure'};

decomp_sexratio_back.Statistic = {'Married women housework';...
    'Married women paid work';'Married women leisure';...
    'Married men housework';'Married men paid work';...
    'Married men leisure'};

ta_married_1990 = struct2array(model.ta1990)'*118;
ta_married_1990 = ta_married_1990(1:6);

ta_married_2010 = struct2array(model.ta2010)'*118;
ta_married_2010 = ta_married_2010(1:6);

ta_married_d1 = struct2array(decomp.tad1)'*118;
ta_married_d1 = ta_married_d1(1:6);

ta_married_d5 = struct2array(decomp.tad5)'*118;
ta_married_d5 = ta_married_d5(1:6);

decomp_sexratio_for.Model_1990 = ta_married_1990;
decomp_sexratio_for.Bargaining = struct2array(decomp.tad51)'*118;
decomp_sexratio_for.Marital_sorting = struct2array(decomp.tad52)'*118;
decomp_sexratio_for.Total_sex_ratio = ta_married_d5;

decomp_sexratio_back.Model_2010 = ta_married_2010;
decomp_sexratio_back.Bargaining = struct2array(decomp.tad11)'*118;
decomp_sexratio_back.Marital_sorting = struct2array(decomp.tad12)'*118;
decomp_sexratio_back.Total_sex_ratio = ta_married_d1;

decomposition_sexratio_fw = struct2table(decomp_sexratio_for);
decomposition_sexratio_bk = struct2table(decomp_sexratio_back);

cd(results_dir);

Sheet = "dsexratio_fw, sigma=" + param.sigma;

writetable(decomposition_sexratio_fw,filename_results,...
    'Sheet',Sheet);

Sheet = "dsexratio_bk, sigma=" + param.sigma;

writetable(decomposition_sexratio_bk,filename_results,...
    'Sheet',Sheet);

%% Quantitative experiment 1: sex ratio goes up to 1.2

cd(main_dir);

% Set the parameters, wages and distributions to 2010 values, but the sex
% ratio to 1.2

param.pe = 0.5;
param.A_g = (1+param.growth_A_g)^20;

param.theta_0 = 1.2;

param.Pf = data.skilldist.Pf2010;
param.Pm = data.skilldist.Pm2010;

wages.m = wages.m2010;
wages.f = wages.f2010;

param.sigma_c = uw_singlefemale(1);
param.sigma_l = uw_singlefemale(2);
param.sigma_g = uw_singlefemale(3);
USfqe1 = uflow_singles(param,wages.f);

param.sigma_c = uw_singlemale(1);
param.sigma_l = uw_singlemale(2);
param.sigma_g = uw_singlemale(3);
USmqe1= uflow_singles(param,wages.m);

param.sigma_c = married_sigma_c;
param.sigma_l = married_sigma_l;
param.sigma_g = 1-param.sigma_c-param.sigma_l;

% Compute equilibrium marriage objects
[theta_S_qe1,Q_qe1,PWf_qe1] = sseq(param,wages,USfqe1,USmqe1,...
    'bargaining','Egalitarian');

MS_qe1 = marital_sorting(param,theta_S_qe1,Q_qe1);
q_experiment.MS.CMqe1 = MS_qe1;

avgPWfqe1 = sum(sum(MS_qe1.*PWf_qe1));

[Cf_qe1,Cm_qe1,Lf_qe1,Lm_qe1,G_qe1,Hf_qe1,Hm_qe1,H_qe1,EQ_qe1,...
    Nf_qe1,Nm_qe1] = per_sol_married(param,wages.f,wages.m,PWf_qe1);
            
q_experiment.taqe1.hfmarried = sum(sum(MS_qe1.*Hf_qe1));
q_experiment.taqe1.nfmarried = sum(sum(MS_qe1.*Nf_qe1));
q_experiment.taqe1.lfmarried = sum(sum(MS_qe1.*Lf_qe1));
q_experiment.taqe1.hmmarried = sum(sum(MS_qe1.*Hm_qe1));
q_experiment.taqe1.nmmarried = sum(sum(MS_qe1.*Nm_qe1));
q_experiment.taqe1.lmmarried = sum(sum(MS_qe1.*Lm_qe1));

% Single female time allocation
param.sigma_c = uw_singlefemale(1);
param.sigma_l = uw_singlefemale(2);
param.sigma_g = uw_singlefemale(3);
[~,lfs,~,hfs,~,nfs] = per_sol_singles(param,wages.f);

q_experiment.taqe1.hfsingle = sum(hfs.*param.Pf);
q_experiment.taqe1.nfsingle = sum(nfs.*param.Pf);
q_experiment.taqe1.lfsingle = sum(lfs.*param.Pf);

% Single men time allocation
param.sigma_c = uw_singlemale(1);
param.sigma_l = uw_singlemale(2);
param.sigma_g = uw_singlemale(3);

[~,lms,~,hms,~,nms] = per_sol_singles(param,wages.m);

q_experiment.taqe1.hmsingle = sum(hms.*param.Pm);
q_experiment.taqe1.nmsingle = sum(nms.*param.Pm);
q_experiment.taqe1.lmsingle = sum(lms.*param.Pm);

% Table for quantitative experiment 1

qe1.Statistic = {'Married women housework';...
    'Married women paid work';'Married women leisure';...
    'Married men housework';'Married men paid work';...
    'Married men leisure';'Single women housework';...
    'Single women paid work';'Single women leisure';...
    'Single men housework';'Single men paid work';...
    'Single men leisure';'Married women consumption';...
    'Married men consumption';'Average wife Pareto weight';...
    'Assortative mating measure'};

qe1.Model2010 = [struct2array(model.ta2010)'*118;...
    sum(Cf_2010.*MS_2010,'all');sum(Cm_2010.*MS_2010,'all');...
    avgPWf2010;msmeasures(model.MS.CM2010)];

qe1.Sex_ratio_12 = [struct2array(q_experiment.taqe1)'*118;...
    sum(Cf_qe1.*MS_qe1,'all');sum(Cm_qe1.*MS_qe1,'all');...
    avgPWfqe1;msmeasures(q_experiment.MS.CMqe1)];

quant_exp_1_results = struct2table(qe1);

cd(results_dir);

Sheet = "quant_exp_1_results, sigma=" + param.sigma;

writetable(quant_exp_1_results,filename_results,'Sheet',Sheet);

%% Quantitative experiment 2: gender wage ratio stays constant

cd(main_dir);

% Set the parameters, male wages and distributions to 2010 values, but the
% gender wage ratio to its 1990 value

param.pe = 0.5;
param.A_g = (1+param.growth_A_g)^20;

param.theta_0 = data.sexratio.sexratio2010;

param.Pf = data.skilldist.Pf2010;
param.Pm = data.skilldist.Pm2010;

wages.m = wages.m2010;
wages.f = wages.m2010*wages.phi1990;

param.sigma_c = uw_singlefemale(1);
param.sigma_l = uw_singlefemale(2);
param.sigma_g = uw_singlefemale(3);
USfqe2 = uflow_singles(param,wages.f);

param.sigma_c = uw_singlemale(1);
param.sigma_l = uw_singlemale(2);
param.sigma_g = uw_singlemale(3);
USmqe2= uflow_singles(param,wages.m);

param.sigma_c = married_sigma_c;
param.sigma_l = married_sigma_l;
param.sigma_g = 1-param.sigma_c-param.sigma_l;

% Compute equilibrium marriage objects
[theta_S_qe2,Q_qe2,PWf_qe2] = sseq(param,wages,USfqe2,USmqe2,...
    'bargaining','Egalitarian');

MS_qe2 = marital_sorting(param,theta_S_qe2,Q_qe2);
q_experiment.MS.CMqe2 = MS_qe2;

avgPWfqe2 = sum(sum(MS_qe2.*PWf_qe2));

[Cf_qe2,Cm_qe2,Lf_qe2,Lm_qe2,G_qe2,Hf_qe2,Hm_qe2,H_qe2,EQ_qe2,...
    Nf_qe2,Nm_qe2] = per_sol_married(param,wages.f,wages.m,PWf_qe2);
            
q_experiment.taqe2.hfmarried = sum(sum(MS_qe2.*Hf_qe2));
q_experiment.taqe2.nfmarried = sum(sum(MS_qe2.*Nf_qe2));
q_experiment.taqe2.lfmarried = sum(sum(MS_qe2.*Lf_qe2));
q_experiment.taqe2.hmmarried = sum(sum(MS_qe2.*Hm_qe2));
q_experiment.taqe2.nmmarried = sum(sum(MS_qe2.*Nm_qe2));
q_experiment.taqe2.lmmarried = sum(sum(MS_qe2.*Lm_qe2));

% Single female time allocation
param.sigma_c = uw_singlefemale(1);
param.sigma_l = uw_singlefemale(2);
param.sigma_g = uw_singlefemale(3);
[~,lfs,~,hfs,~,nfs] = per_sol_singles(param,wages.f);

q_experiment.taqe2.hfsingle = sum(hfs.*param.Pf);
q_experiment.taqe2.nfsingle = sum(nfs.*param.Pf);
q_experiment.taqe2.lfsingle = sum(lfs.*param.Pf);

% Single men time allocation
param.sigma_c = uw_singlemale(1);
param.sigma_l = uw_singlemale(2);
param.sigma_g = uw_singlemale(3);

[~,lms,~,hms,~,nms] = per_sol_singles(param,wages.m);

q_experiment.taqe2.hmsingle = sum(hms.*param.Pm);
q_experiment.taqe2.nmsingle = sum(nms.*param.Pm);
q_experiment.taqe2.lmsingle = sum(lms.*param.Pm);

% Table for quantitative experiment 2

qe2.Statistic = {'Married women housework';...
    'Married women paid work';'Married women leisure';...
    'Married men housework';'Married men paid work';...
    'Married men leisure';'Single women housework';...
    'Single women paid work';'Single women leisure';...
    'Single men housework';'Single men paid work';...
    'Single men leisure';'Married women consumption';...
    'Married men consumption';'Average wife Pareto weight';...
    'Assortative mating measure'};

qe2.Data2010 = [struct2array(data.ta2010)'*118;...
    nan ; nan ; nan ; msmeasures(data.MS.CM2010)];

qe2.Model2010 = [struct2array(model.ta2010)'*118;...
    sum(Cf_2010.*MS_2010,'all');sum(Cm_2010.*MS_2010,'all');...
    avgPWf2010;msmeasures(model.MS.CM2010)];

qe2.Constant_gender_wage_ratio = [struct2array(q_experiment.taqe2)'*118;...
    sum(Cf_qe2.*MS_cf1,'all');sum(Cm_qe2.*MS_qe2,'all');...
    avgPWfqe2;msmeasures(q_experiment.MS.CMqe2)];

quant_exp_2_results = struct2table(qe2);

cd(results_dir);

Sheet = "quant_exp_2_results, sigma=" + param.sigma;

writetable(quant_exp_2_results,filename_results,'Sheet',Sheet);

%% Quantitative experiment 3: no gender wage ratio 

cd(main_dir);

% Set the parameters, male wages and distributions to 2010 values, but the
% gender wage ratio to its 1990 value

param.pe = 0.5;
param.A_g = (1+param.growth_A_g)^20;

param.theta_0 = data.sexratio.sexratio2010;

param.Pf = data.skilldist.Pf2010;
param.Pm = data.skilldist.Pm2010;

wages.m = wages.m2010;
wages.f = wages.m;

param.sigma_c = uw_singlefemale(1);
param.sigma_l = uw_singlefemale(2);
param.sigma_g = uw_singlefemale(3);
USfqe3 = uflow_singles(param,wages.f);

param.sigma_c = uw_singlemale(1);
param.sigma_l = uw_singlemale(2);
param.sigma_g = uw_singlemale(3);
USmqe3= uflow_singles(param,wages.m);

param.sigma_c = married_sigma_c;
param.sigma_l = married_sigma_l;
param.sigma_g = 1-param.sigma_c-param.sigma_l;

% Compute equilibrium marriage objects
[theta_S_qe3,Q_qe3,PWf_qe3] = sseq(param,wages,USfqe3,USmqe3,...
    'bargaining','Egalitarian');

MS_qe3 = marital_sorting(param,theta_S_qe3,Q_qe3);
q_experiment.MS.CMqe3 = MS_qe3;

avgPWfqe3 = sum(sum(MS_qe3.*PWf_qe3));

[Cf_qe3,Cm_qe3,Lf_qe3,Lm_qe3,G_qe3,Hf_qe3,Hm_qe3,H_qe3,EQ_qe3,...
    Nf_qe3,Nm_qe3] = per_sol_married(param,wages.f,wages.m,PWf_qe3);
            
q_experiment.taqe3.hfmarried = sum(sum(MS_qe3.*Hf_qe3));
q_experiment.taqe3.nfmarried = sum(sum(MS_qe3.*Nf_qe3));
q_experiment.taqe3.lfmarried = sum(sum(MS_qe3.*Lf_qe3));
q_experiment.taqe3.hmmarried = sum(sum(MS_qe3.*Hm_qe3));
q_experiment.taqe3.nmmarried = sum(sum(MS_qe3.*Nm_qe3));
q_experiment.taqe3.lmmarried = sum(sum(MS_qe3.*Lm_qe3));

% Single female time allocation
param.sigma_c = uw_singlefemale(1);
param.sigma_l = uw_singlefemale(2);
param.sigma_g = uw_singlefemale(3);
[~,lfs,~,hfs,~,nfs] = per_sol_singles(param,wages.f);

q_experiment.taqe3.hfsingle = sum(hfs.*param.Pf);
q_experiment.taqe3.nfsingle = sum(nfs.*param.Pf);
q_experiment.taqe3.lfsingle = sum(lfs.*param.Pf);

% Single men time allocation
param.sigma_c = uw_singlemale(1);
param.sigma_l = uw_singlemale(2);
param.sigma_g = uw_singlemale(3);

[~,lms,~,hms,~,nms] = per_sol_singles(param,wages.m);

q_experiment.taqe3.hmsingle = sum(hms.*param.Pm);
q_experiment.taqe3.nmsingle = sum(nms.*param.Pm);
q_experiment.taqe3.lmsingle = sum(lms.*param.Pm);

% Table for quantitative experiment 3

qe3.Statistic = {'Married women housework';...
    'Married women paid work';'Married women leisure';...
    'Married men housework';'Married men paid work';...
    'Married men leisure';'Single women housework';...
    'Single women paid work';'Single women leisure';...
    'Single men housework';'Single men paid work';...
    'Single men leisure';'Married women consumption';...
    'Married men consumption';'Average wife Pareto weight';...
    'Assortative mating measure'};

qe3.Model2010 = [struct2array(model.ta2010)'*118;...
    sum(Cf_2010.*MS_2010,'all');sum(Cm_2010.*MS_2010,'all');...
    avgPWf2010;msmeasures(model.MS.CM2010)];

qe3.No_gender_wage_gap_2010 = [struct2array(q_experiment.taqe3)'*118;...
    sum(Cf_qe3.*MS_cf1,'all');sum(Cm_qe3.*MS_qe3,'all');...
    avgPWfqe3;msmeasures(q_experiment.MS.CMqe3)];

quant_exp_3_results = struct2table(qe3);

cd(results_dir);

Sheet = "quant_exp_3_results, sigma=" + param.sigma;

writetable(quant_exp_3_results,filename_results,'Sheet',Sheet);

% %% Quantitative exercise 1
% 
% param.theta_0 = data.sexratio_90;
% 
% [theta_S_sim1,Q_sim1,PWf_sim1] = sseq(param,wages,USf10,USm10,...
%     'bargaining','Egalitarian',...
%     'quietly','true');
% 
% MS_sim1 = marital_sorting(param,theta_S_sim1,Q_sim1);
% sim1.MS.MSsim1 = MS_sim1;
% avgPWfsim1 = sum(sum(MS_sim1.*PWf_sim1));
% 
% randmat_sim1 = sum(MS_sim1,2)*sum(MS_sim1,1);
% msdelta_sim1 = trace(MS_sim1)/trace(randmat_sim1);
% 
% [Cf_sim1,Cm_sim1,Lf_sim1,Lm_sim1,G_sim1,Hf_sim1,Hm_sim1,H_sim1,EQ_sim1,...
%     Nf_sim1,Nm_sim1] = per_sol_married(param,wages.f,wages.m,PWf_sim1);
%          
% sim1.ta.hfmarried = sum(sum(MS_sim1.*Hf_sim1));
% sim1.ta.nfmarried = sum(sum(MS_sim1.*Nf_sim1));
% sim1.ta.lfmarried = sum(sum(MS_sim1.*Lf_sim1));
% sim1.ta.hmmarried = sum(sum(MS_sim1.*Hm_sim1));
% sim1.ta.nmmarried = sum(sum(MS_sim1.*Nm_sim1));
% sim1.ta.lmmarried = sum(sum(MS_sim1.*Lm_sim1));
% 
% %% Quantitative exercise 2
% 
% param.theta_0 = 1.2;
% 
% [theta_S_sim2,Q_sim2,PWf_sim2] = sseq(param,wages,USf10,USm10,...
%     'bargaining','Egalitarian',...
%     'quietly','true');
% 
% MS_sim2 = marital_sorting(param,theta_S_sim2,Q_sim2);
% sim2.MS.MSsim2 = MS_sim2;
% avgPWfsim2 = sum(sum(MS_sim2.*PWf_sim2));
% 
% randmat_sim2 = sum(MS_sim2,2)*sum(MS_sim2,1);
% msdelta_sim2 = trace(MS_sim2)/trace(randmat_sim2);
% 
% [Cf_sim2,Cm_sim2,Lf_sim2,Lm_sim2,G_sim2,Hf_sim2,Hm_sim2,H_sim2,EQ_sim2,...
%     Nf_sim2,Nm_sim2] = per_sol_married(param,wages.f,wages.m,PWf_sim2);
%                   
% sim2.ta.hfmarried = sum(sum(MS_sim2.*Hf_sim2));
% sim2.ta.nfmarried = sum(sum(MS_sim2.*Nf_sim2));
% sim2.ta.lfmarried = sum(sum(MS_sim2.*Lf_sim2));
% sim2.ta.hmmarried = sum(sum(MS_sim2.*Hm_sim2));
% sim2.ta.nmmarried = sum(sum(MS_sim2.*Nm_sim2));
% sim2.ta.lmmarried = sum(sum(MS_sim2.*Lm_sim2));
% 
% %% Quantitative exercise 3
% 
% param.theta_0 = data.sexratio_10;
% wages.f = wages.m;
% 
% [theta_S_sim3,Q_sim3,PWf_sim3] = sseq(param,wages,USf10,USm10,...
%     'bargaining','Egalitarian',...
%     'quietly','true');
% 
% MS_sim3 = marital_sorting(param,theta_S_sim3,Q_sim3);
% sim3.MS.MSsim3 = MS_sim3;
% avgPWfsim3 = sum(sum(MS_sim3.*PWf_sim3));
% 
% randmat_sim3 = sum(MS_sim3,2)*sum(MS_sim3,1);
% msdelta_sim3 = trace(MS_sim3)/trace(randmat_sim3);
% 
% [Cf_sim3,Cm_sim3,Lf_sim3,Lm_sim3,G_sim3,Hf_sim3,Hm_sim3,H_sim3,EQ_sim3,...
%     Nf_sim3,Nm_sim3] = per_sol_married(param,wages.f,wages.m,PWf_sim3);
%                   
% sim3.ta.hfmarried = sum(sum(MS_sim3.*Hf_sim3));
% sim3.ta.nfmarried = sum(sum(MS_sim3.*Nf_sim3));
% sim3.ta.lfmarried = sum(sum(MS_sim3.*Lf_sim3));
% sim3.ta.hmmarried = sum(sum(MS_sim3.*Hm_sim3));
% sim3.ta.nmmarried = sum(sum(MS_sim3.*Nm_sim3));
% sim3.ta.lmmarried = sum(sum(MS_sim3.*Lm_sim3));
% 
% %% Export the results
% 
% % Set up path and filename
% 
% results_path = 'G:\My Drive\China project\Mannheim seminar';
% cd(results_path);
% filename_results = 'Results_Mannheim.xlsx';
% 
% % Time allocation results
% 
% ta_results.Statistic = {'Wives housework';'Wives paid work';...
%     'Wives leisure';'Husbands housework';'Husbands paid work';...
%     'Husbands leisure';'Single women housework';...
%     'Single women paid work';'Single women leisure';...
%     'Single men housework';'Single men paid work';...
%     'Single men leisure'};
% 
% ta_results.Data90 = struct2array(data.ta90)'*118;
% ta_results.Model90 = struct2array(model.ta90)'*118;
% ta_results.Data10 = struct2array(data.ta10)'*118;
% ta_results.Model10 = struct2array(model.ta10)'*118;
% 
% time_allocation_results = struct2table(ta_results);
% 
% writetable(time_allocation_results,filename_results,...
%     'Sheet','Time allocation results');
% 
% % Time allocation within wife/husband type cells
% writematrix(Nf_90,filename_results,'Sheet','Model detailed LFP 1990');
% writematrix(Nf_10,filename_results,'Sheet','Model detailed LFP 2010');
% % Nf90
% % Nf10
% 
% % Marital sorting results
% 
% % Obtain the contingency matrices 
% data.MS.sortcontmat90 = data.MS.MS90/sum(data.MS.MS90,'all');
% data.MS.sortcontmat10 = data.MS.MS10/sum(data.MS.MS10,'all');
% 
% model.MS.sortcontmat90 = model.MS.MS90/sum(model.MS.MS90,'all');
% model.MS.sortcontmat10 = model.MS.MS10/sum(model.MS.MS10,'all');
% 
% % Obtain the random matrices
% data.MS.randmat90 = ...
%     sum(data.MS.sortcontmat90,2)*sum(data.MS.sortcontmat90,1);
% data.MS.randmat10 = ...
%     sum(data.MS.sortcontmat10,2)*sum(data.MS.sortcontmat10,1);
% 
% model.MS.randmat90 = ...
%     sum(model.MS.MS90,2)*sum(model.MS.MS90,1);
% model.MS.randmat10 = ...
%     sum(model.MS.MS10,2)*sum(model.MS.MS10,1);
% 
% % Obtain the marital sorting measures
% data.MS.msdelta90 = trace(data.MS.sortcontmat90)/trace(data.MS.randmat90);
% data.MS.msdelta10 = trace(data.MS.sortcontmat10)/trace(data.MS.randmat10);
% 
% model.MS.msdelta90 = ...
%     trace(model.MS.MS90)/trace(model.MS.randmat90);
% model.MS.msdelta10 = ...
%     trace(model.MS.MS10)/trace(model.MS.randmat10);
% 
% % Export everything
% 
% writematrix(data.MS.sortcontmat90,filename_results,...
%     'Sheet','Data MS 1990');
% writematrix(data.MS.randmat90,filename_results,...
%     'Sheet','Data random MS 1990');
% writematrix(data.MS.sortcontmat10,filename_results,...
%     'Sheet','Data MS 2010');
% writematrix(data.MS.randmat10,filename_results,...
%     'Sheet','Data random MS 2010');
% 
% writematrix(model.MS.sortcontmat90,filename_results,...
%     'Sheet','Model MS 1990');
% writematrix(model.MS.randmat90,filename_results,...
%     'Sheet','Model random MS 1990');
% writematrix(model.MS.sortcontmat10,filename_results,...
%     'Sheet','Model MS 2010');
% writematrix(model.MS.randmat10,filename_results,...
%     'Sheet','Model random MS 2010');
% 
% % Simulations 
% 
% sim_results.Statistic = {'Wives housework';'Wives paid work';...
%     'Wives leisure';'Husbands housework';'Husbands paid work';...
%     'Husbands leisure';'Assortative mating measure';...
%     'Average Pareto weight wives'};
% 
% modelta10 = struct2array(model.ta10)'*118;
% 
% sim_results.Model2010 = [modelta10(1:6);model.MS.msdelta10;avgPWf10];
% 
% ta_sim1 = struct2array(sim1.ta)'*118;
% ta_sim2 = struct2array(sim2.ta)'*118;
% ta_sim3 = struct2array(sim3.ta)'*118;
% 
% sim_results.Sim1 = [ta_sim1;msdelta_sim1;avgPWfsim1];
% sim_results.Sim2 = [ta_sim2;msdelta_sim2;avgPWfsim2];
% sim_results.Sim3 = [ta_sim3;msdelta_sim3;avgPWfsim3];
% 
% sim_results_table = struct2table(sim_results);
% 
% writetable(sim_results_table,filename_results,...
%     'Sheet','Simulation results');
% 
% % Means of match quality draws
% 
% writematrix(param.MU,filename_results,...
%     'Sheet','Means MQD');
% 
% % Consumption changes
% 
% cons.Skill_level = {'Low';'Medium';'High'};
% 
% cons.Model10 = Cf_10*param.Pf;
% cons.Sim1 = Cf_sim1*param.Pf;
% cons.Sim2 = Cf_sim2*param.Pf;
% cons.Sim3 = Cf_sim3*param.Pf;
% 
% cons_table = struct2table(cons);
% 
% writetable(cons_table,filename_results,...
%     'Sheet','Consumption by skill level');
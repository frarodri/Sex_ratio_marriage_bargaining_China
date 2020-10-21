%{
Calibration of the whole model: this code computes the set of parameters
that minimize the distance between the data targets and the model
counterparts. There are two main groups of targets: marital sorting and
time allocation. The parameters that discipline the time allocation moments
are those in the utility function of married people. The parameters
controlling the marital sorting are the means of the match quality draws
among entrants of different types.
%} 

% home: G:\My Drive\China project\Model\Sub-markets  
% uni: G:\Mi unidad\China project\Model\Sub-markets

clear ; close all; clc
cd('G:\My Drive\China project\Model\Sub-markets')
%%
%%%%%%%%%%%%%%%%%
% Read the data %
%%%%%%%%%%%%%%%%%

filename = 'Data_for_model.xlsx';

% Sex ratio of young population
data.theta_0 = readmatrix(filename,'Sheet','Sex ratio','Range','B2:B3');

% Education distribution
data.Pf = ...
    readmatrix(filename,'Sheet','Education females','Range','B2:F3')';

data.Pm = ...
    readmatrix(filename,'Sheet','Education males','Range','B2:F3')';

% Marital sorting
data.contmat_1990 = ...
    readmatrix(filename,'Sheet','Marital Sorting','Range','B3:F7');
data.contmat_2010 = ...
    readmatrix(filename,'Sheet','Marital Sorting','Range','H3:L7');

for i = {'contmat_1990','contmat_2010'}
   data.(i{1}) = data.(i{1})/sum(sum(data.(i{1}))); 
end

% Wages
data.gender_wage_ratio = ...
    readmatrix(filename,'Sheet','Gender wage ratio','Range','B2:B3');
data.skill_premium_2010 = ...
    readmatrix(filename,'Sheet','Skill premium','Range','B2:B6');
data.wage_growth = ...
    readmatrix(filename,'Sheet','Wage growth','Range','A2:A2');

data.wagem_1990 = ones(size(data.skill_premium_2010));
data.wagef_1990 = data.wagem_1990*data.gender_wage_ratio(1);

g = data.wage_growth/(data.Pm(:,2)'*data.skill_premium_2010);
data.wagem_2010 = g*data.skill_premium_2010;
data.wagef_2010 = data.wagem_2010*data.gender_wage_ratio(2);

% Time allocation
data.time_allocation = ...
    readmatrix(filename,'Sheet','Time allocation','Range','D2:F9');

% Prices for home equipment and home production productivity
%data.pe_growth = 
%data.homeprodgrowth = 
data.home_equip_share = 0.05; % U.S. number in Knowles (2014)

% Life expectancy
data.life_expectancy = [69.293;75.236]; % World bank

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Parameters chosen externally that are constant %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

param.beta = 0.96; % Standard
param.sigma = 1.5; % Attanasio et. al (2008)
param.alpha_x = 1/2; % Makes sense for a matching function

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Parameters that change in time for 1990 %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

lifexpect90 = data.life_expectancy(1);
param.delta = 1/(lifexpect90-20);

param.pe = 1;
param.A_g = 1;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Elasticity of home production with respect to home equipment %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
param.alpha_g = 1-data.home_equip_share;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Weights in the utility functions of singles %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Define the inputs for the minimization routine
x0 = 1/3*[1 1 1];
A = [];
b = []; 
Aeq = [1 1 1];
beq = 1;
lb = [0 0 0];
ub = [];

options = optimoptions('fmincon','Display','off');

% Extract the time allocation data 
hsingf90 = data.time_allocation(3,1);
nsingf90 = data.time_allocation(3,2);
lsingf90 = data.time_allocation(3,3);

hsingm90 = data.time_allocation(1,1);
nsingm90 = data.time_allocation(1,2);
lsingm90 = data.time_allocation(1,3);

wagef90 = data.gender_wage_ratio(1);
wagem90 = 1;

% Single women
% Define the targets 
targ.h = hsingf90;
targ.n = nsingf90;
targ.l = lsingf90;

% Set the wage
wage = wagef90;

% Minimize the loss function
uweightsingf = fmincon(@(x)timealloclossfn_s(param,targ,wage,x),x0,A,b,...
                                                Aeq,beq,lb,ub,[],options);
                                            
uweightsing.c_f = uweightsingf(1); 
uweightsing.l_f = uweightsingf(2); 
uweightsing.g_f = uweightsingf(3); 
                                            
% Single men
% Define the targets
targ.h = hsingm90;
targ.n = nsingm90;
targ.l = lsingm90;

% Set the wage
wage = wagem90;

% Minimize the loss function
uweightsingm = fmincon(@(x)timealloclossfn_s(param,targ,wage,x),x0,A,b,...
                                                Aeq,beq,lb,ub,[],options);

uweightsing.c_m = uweightsingm(1); 
uweightsing.l_m = uweightsingm(2); 
uweightsing.g_m = uweightsingm(3);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Parameters of the housework time aggregator for married couples %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Define gender housework and wage ratios for each year
housewrkrat90 = data.time_allocation(4,1)/data.time_allocation(2,1);
housewrkrat10 = data.time_allocation(8,1)/data.time_allocation(6,1);

gwagerat90 = data.gender_wage_ratio(1);
gwagerat10 = data.gender_wage_ratio(2);

% param.eta = (log(1/gwagerat10)-log(1/gwagerat90))/...
%             (log(housewrkrat10)-log(housewrkrat90));

param.eta = 0.33;

param.eta_f = (housewrkrat90^param.eta)*gwagerat90/...
                (1+(housewrkrat90^param.eta)*gwagerat90);
            
%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Now, the big loop to find the utility weights for married people, the %
% means of the match quality draws and the search frictions in the      %
% same-type marriage markets                                            %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% For now, we set (later we need to make this automatic, maybe bisection
% method...)

param.A_x = 1;
param.sigma_c = 0.343;
param.sigma_l = 0.618; % For sigma=1.5
param.sigma_g = 0.039;
% param.sigma_c = 0.375;
% param.sigma_l = 0.535; % For sigma=1
% param.sigma_g = 0.09;

% param.singcostf = 1.185; % NEW PARAMETER, COST OF BEING SINGLE FOR WOMEN
param.singcostf = 1.29;

% Compute the Steady State equilibrium and the means for the match quality
% draws

% Set sex ratio, distributions over types and wages to 1990 level
param.theta_0 = data.theta_0(1);
param.Pf = data.Pf(:,1);
param.Pm = data.Pm(:,1);

wages.m = data.wagem_1990;
wages.f = data.wagef_1990;

targ.MS = data.contmat_1990;
targ.MS(5,1) = 10^(-8);
targ.MS = targ.MS/sum(sum(targ.MS));

% [theta_S,Q,PWf,MU] = sseq_meandfind(param,wages,targ);

param.MU = zeros(5,5);
[theta_S_90,Q_90,PWf_90] = sseq(param,wages);

types = size(Q_90,1); 
x0 = zeros(types*(types-1),1);

options = optimoptions('fminunc','Display','off');
mus = fminunc(@(MU_vec)l_maritalsort(param,targ,theta_S_90,Q_90,MU_vec),...
    x0,options);

MU_aux = vec2mat(mus,types);
MU_l = [zeros(1,types);tril(MU_aux)];
MU_u = [triu(MU_aux,1);zeros(1,types)];
MU = MU_l+MU_u;

param.MU = MU;
MS_90 = marital_sorting(param,theta_S_90,Q_90);
avgPWf_90 = sum(sum(PWf_90.*MS_90)); 

[Cf_90,Cm_90,Lf_90,Lm_90,G_90,Hf_90,Hm_90,H_90,EQ_90,Nf_90,Nm_90] = ...
                      per_sol_married(param,wages.f,wages.m,PWf_90);
                 
lfmarried90 = sum(sum(MS_90.*Lf_90));
lmmarried90 = sum(sum(MS_90.*Lm_90));
hfmarried90 = sum(sum(MS_90.*Hf_90));
hmmarried90 = sum(sum(MS_90.*Hm_90));
nfmarried90 = sum(sum(MS_90.*Nf_90));
nmmarried90 = sum(sum(MS_90.*Nm_90));

%%
%%%%%%%%%%%%%%%%%%%%%%
% See how 2011 fares %
%%%%%%%%%%%%%%%%%%%%%%

% Set sex ratio, distributions over types and wages to 2011 level
param.theta_0 = data.theta_0(2);
param.Pf = data.Pf(:,2);
param.Pm = data.Pm(:,2);

wages.m = data.skill_premium_2010;
wages.f = wages.m*data.gender_wage_ratio(2);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Parameters that change in time for 1990 %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

lifexpect2010 = data.life_expectancy(2);
param.delta = 1/(lifexpect2010-20);

param.pe = 1/3;
param.A_g = g;


[theta_S_2010,Q_2010,PWf_2010] = sseq(param,wages);

MS_2010 = marital_sorting(param,theta_S_2010,Q_2010);
avgPWf2010 = sum(sum(MS_2010.*PWf_2010));

[Cf_2010,Cm_2010,Lf_2010,Lm_2010,G_2010,Hf_2010,Hm_2010,H_2010,EQ_2010,...
    Nf_2010,Nm_2010] = per_sol_married(param,wages.f,wages.m,PWf_2010);
                  
lfmarried2010 = sum(sum(MS_2010.*Lf_2010));
lmmarried2010 = sum(sum(MS_2010.*Lm_2010));
hfmarried2010 = sum(sum(MS_2010.*Hf_2010));
hmmarried2010 = sum(sum(MS_2010.*Hm_2010));
nfmarried2010 = sum(sum(MS_2010.*Nf_2010));
nmmarried2010 = sum(sum(MS_2010.*Nm_2010));

%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Some preliminary simulations %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%
% What would happen in 2010 had the sex ratio stayed at the 1990 level

param.theta_0 = data.theta_0(1);

[theta_S_sim1,Q_sim1,PWf_sim1] = sseq(param,wages);

MS_sim1 = marital_sorting(param,theta_S_sim1,Q_sim1);
avgPWfsim1 = sum(sum(MS_sim1.*PWf_sim1));

[Cf_sim1,Cm_sim1,Lf_sim1,Lm_sim1,G_sim1,Hf_sim1,Hm_sim1,H_sim1,EQ_sim1,...
    Nf_sim1,Nm_sim1] = per_sol_married(param,wages.f,wages.m,PWf_sim1);
                  
lfmarriedsim1 = sum(sum(MS_sim1.*Lf_sim1));
lmmarriedsim1 = sum(sum(MS_sim1.*Lm_sim1));
hfmarriedsim1 = sum(sum(MS_sim1.*Hf_sim1));
hmmarriedsim1 = sum(sum(MS_sim1.*Hm_sim1));
nfmarriedsim1 = sum(sum(MS_sim1.*Nf_sim1));
nmmarriedsim1 = sum(sum(MS_sim1.*Nm_sim1));

%%
% What would happen in 2010 had the sex ratio been already at 1.2

param.theta_0 = 1.2;

[theta_S_sim2,Q_sim2,PWf_sim2] = sseq(param,wages);

MS_sim2 = marital_sorting(param,theta_S_sim2,Q_sim2);
avgPWfsim2 = sum(sum(MS_sim2.*PWf_sim2));

[Cf_sim2,Cm_sim2,Lf_sim2,Lm_sim2,G_sim2,Hf_sim2,Hm_sim2,H_sim2,EQ_sim2,...
    Nf_sim2,Nm_sim2] = per_sol_married(param,wages.f,wages.m,PWf_sim2);
                  
lfmarriedsim2 = sum(sum(MS_sim2.*Lf_sim2));
lmmarriedsim2 = sum(sum(MS_sim2.*Lm_sim2));
hfmarriedsim2 = sum(sum(MS_sim2.*Hf_sim2));
hmmarriedsim2 = sum(sum(MS_sim2.*Hm_sim2));
nfmarriedsim2 = sum(sum(MS_sim2.*Nf_sim2));
nmmarriedsim2 = sum(sum(MS_sim2.*Nm_sim2));

%%
% What would happen in 2010 with a higher sex ratio AND reduced gender wage
% gap (gender wage ratio 0.85)

wages.f = wages.m;

[theta_S_sim3,Q_sim3,PWf_sim3] = sseq(param,wages);

MS_sim3 = marital_sorting(param,theta_S_sim3,Q_sim3);
avgPWfsim3 = sum(sum(MS_sim3.*PWf_sim3));

[Cf_sim3,Cm_sim3,Lf_sim3,Lm_sim3,G_sim3,Hf_sim3,Hm_sim3,H_sim3,EQ_sim3,...
    Nf_sim3,Nm_sim3] = per_sol_married(param,wages.f,wages.m,PWf_sim3);
                  
lfmarriedsim3 = sum(sum(MS_sim3.*Lf_sim3));
lmmarriedsim3 = sum(sum(MS_sim3.*Lm_sim3));
hfmarriedsim3 = sum(sum(MS_sim3.*Hf_sim3));
hmmarriedsim3 = sum(sum(MS_sim3.*Hm_sim3));
nfmarriedsim3 = sum(sum(MS_sim3.*Nf_sim3));
nmmarriedsim3 = sum(sum(MS_sim3.*Nm_sim3));

%%
% Tables with results

Statistic = {'Wives housework time' ; 'Wives paid work time' ; ...
    'Wives leisure time' ; 'Husbands housework time' ; ...
    'Husbands paid work time' ; 'Husbands leisure time'};

ta_data_1990 = [data.time_allocation(4,:)' ; data.time_allocation(2,:)']...
    *118;

ta_model_1990 = [hfmarried90 ; nfmarried90 ; lfmarried90 ; ...
    hmmarried90 ; nmmarried90 ; lmmarried90]*118;

varnames_ta = {'Statistic' , 'Data' , 'Model'}; 

T_time_allocation_results_1990 = ...
    table(Statistic,ta_data_1990,ta_model_1990,...
    'VariableNames',varnames_ta);

  


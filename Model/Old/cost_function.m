function [cost,model_mom] = ...
    cost_function(data_mom,param,expect,wages,param_to_calib)

param.lambda = param_to_calib(1);
param.cbar = param_to_calib(2);
param.kappa = param_to_calib(3);
param.A = param_to_calib(4);

param.MU = param.kappa*eye(length(expect.P_f));

[Q, ~, ~, ~, MS] = marriage_pree(param,expect,wages);

[pi_f,pi_m] = match_prob(param,expect.theta_S);

mr_f = pi_f*expect.P_f'*(1-normcdf(Q,param.MU))*expect.P_m;
mr_m = pi_m*expect.P_m'*(1-normcdf(Q,param.MU))'*expect.P_f;
                
model_mom.exp_single_f = 1/mr_f;
model_mom.exp_single_m = 1/mr_m;

[model_mom.delta, ~, ~] = sort_cont_mat(MS);

Lf = married_fls(param,wages);
model_mom.lf_bar = sum(sum(Lf.*MS));

cost = (1/data_mom.delta)^2*(model_mom.delta-data_mom.delta)^2+...
    (1/data_mom.exp_single_f)^2*...
    (model_mom.exp_single_f-data_mom.exp_single_f)^2+...
    (1/data_mom.exp_single_m)^2*...
    (model_mom.exp_single_m-data_mom.exp_single_m)^2+...
    (1/data_mom.lf_bar)^2*(model_mom.lf_bar-data_mom.lf_bar)^2;




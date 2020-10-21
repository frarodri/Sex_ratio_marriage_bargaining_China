function obj = obj_func(param,expect,wages,targets)

[Q, ~, ~, ~, MS] = ...
    marriage_pree(param,expect,wages);

[delta, ~, ~] = sort_cont_mat(MS);

[pi_f,pi_m] = match_prob(param,expect.theta_S);

mr_f = pi_f*expect.Z_f'*(1-normcdf(Q,param.MU))*expect.Z_m;
mr_m = pi_m*expect.Z_f'*(1-normcdf(Q,param.MU))*expect.Z_m;
                
exp_single_f = 1/mr_f;
exp_single_m = 1/mr_m;

obj = (1/targets.delta)*(targets.delta-delta)^2+...
      (1/targets.exp_single_f)*(targets.exp_single_f-exp_single_f)^2+...
      (1/targets.exp_single_m)*(targets.exp_single_m-exp_single_m)^2;


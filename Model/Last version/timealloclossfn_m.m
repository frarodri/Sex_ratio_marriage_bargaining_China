function L = timealloclossfn_m(param,wages,targ,uweightsing,x)

param.sigma_c = x(1);
param.sigma_l = x(2);
param.sigma_g = x(3);

[theta_S,Q,PWf] = sseeq_meanfind(param,wages,targ,uweightsing);

MS = marital_sorting(param,theta_S,Q);

[~,~,Lf,~,~,Hf,~,~,~,Nf,~] = ...
                      per_sol_married(param,wages_f,wages_m,PWf);
                  
L = (sum(sum(MS.*Lf))-targ.lf)^2+(sum(sum(MS.*Hf))-targ.hf)^2+...
    (sum(sum(MS.*Nf))-targ.nf)^2;
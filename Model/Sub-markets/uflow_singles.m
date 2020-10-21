%{
Computes the flow value of single agents and records it in two vectors, one 
for each sex and each having as many entries as types there are for that 
sex
%}

function [USf,USm] = uflow_singles(param,wages)

if iscolumn(wages.f) && iscolumn(wages.m)
       
      % Manual assignment of parameters BAD, fix it
%       param.sigma_c = 0.4035; 
%       param.sigma_l = 0.5306;
%       param.sigma_g = 0.0659;
      
%       param.sigma_c = 0.3526; 
%       param.sigma_l = 0.6229;
%       param.sigma_g = 0.0245;
      
      [Cf,Lf,Gf,~,~,~] = per_sol_singles(param,wages.f);
      USf = ind_u(param,Cf,Lf,Gf)-param.singcostf;
      
%       param.sigma_c = 0.4022; 
%       param.sigma_l = 0.5829;
%       param.sigma_g = 0.0149;
      
%       param.sigma_c = 0.3625; 
%       param.sigma_l = 0.6342;
%       param.sigma_g = 0.0023;
      
      [Cm,Lm,Gm,~,~,~] = per_sol_singles(param,wages.m);      
      USm = ind_u(param,Cm,Lm,Gm);  

else
    
    error('Wage vectors must be columns')
    
end

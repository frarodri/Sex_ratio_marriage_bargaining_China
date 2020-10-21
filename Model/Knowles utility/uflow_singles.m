%{
Computes the flow value of single agents and records it in two vectors, one 
for each sex and each having as many entries as types there are for that 
sex
%}

function [USf,USm] = uflow_singles(param,wages)

if iscolumn(wages.f) && iscolumn(wages.m)
    
      [Cf,Lf,Gf,~,~,~] = per_sol_singles(param,wages.f);
      [Cm,Lm,Gm,~,~,~] = per_sol_singles(param,wages.m);
      
      USf = ind_u(param,Cf,Lf,Gf);
      USm = ind_u(param,Cm,Lm,Gm);  

else
    
    error('Wage vectors must be columns')
    
end

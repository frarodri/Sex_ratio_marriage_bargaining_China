%{
Computes the flow value of single agents and records it in two vectors, one 
for each sex and each having as many entries as types there are for that 
sex
%}

function US = uflow_singles(param,wage_vec)

if iscolumn(wage_vec) && all(wage_vec>0)
      
      [Cf,Lf,Gf,~,~,~] = per_sol_singles(param,wage_vec);
      US = ind_u(param,Cf,Lf,Gf);
      
else
    
    message = ['Vector of wages must be a column with strictly positive '...
        'entries'];
    error(message)
    
    
end

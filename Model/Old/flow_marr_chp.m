%{
Computes the part of the flow value of marriage that corresponds to 
utility of market goods consumption and household production and records 
it in a matrix where rows are the wife's type and columns are the 
husband's type
%}

function FM_CHP = flow_marr_chp(param,wages)

if iscolumn(wages.OMEGA_f) && iscolumn(wages.OMEGA_m)
    
    szOMEGA_f = size(wages.OMEGA_f,1);
    szOMEGA_m = size(wages.OMEGA_m,1);

    Lf = married_fls(param,wages);
    
    C = ones(szOMEGA_f,szOMEGA_m).*wages.OMEGA_m'+Lf.*wages.OMEGA_f;

    N = 1-Lf;

    if param.lambda==1
        
        FM_CHP = log(C-param.gamma*param.cbar);
        
    else
        
        FM_CHP = param.lambda*log(C-param.gamma*param.cbar)...
            +(1-param.lambda)*log(N);
        
    end

else
    
    error('Wage vectors must be columns')
    
end
    
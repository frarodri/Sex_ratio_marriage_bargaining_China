%{
Computes the labor supply of married females and records it in a matrix 
where rows are wife's type and columns husband's type
%}

function Lf = married_fls(param,wages)

if iscolumn(wages.OMEGA_f) && iscolumn(wages.OMEGA_m)
    
    szOMEGA_f = size(wages.OMEGA_f,1);
    szOMEGA_m = size(wages.OMEGA_m,1);

    PHI = wages.OMEGA_f*(1./wages.OMEGA_m)';
    
    Lf = param.lambda-(1-param.lambda)./PHI+...
        (1-param.lambda)*param.gamma*param.cbar*...
        repmat(1./wages.OMEGA_f,1,szOMEGA_m);
        
    Lf = max(zeros(szOMEGA_f,szOMEGA_m),Lf);
    Lf = min(ones(szOMEGA_f,szOMEGA_m),Lf);
    
else
    
    error('Wage vectors must be columns')
    
end
function [VS, ID] = max_single(VS_f,VS_m)

%{
Compute the maximum value of being single and record it in a matrix where 
rows are the wife's type and columns are the husband's type
%}

if iscolumn(VS_f) && iscolumn(VS_m)
    
    szVS_f = size(VS_f,1);
    szVS_m = size(VS_m,1);
    
    VS = zeros(szVS_f,szVS_m);
    ID = zeros(szVS_f,szVS_m);
    eps = exp(-10);
    
    for i = 1:szVS_f
        for j = 1:szVS_m
            
            VS(i,j) = max(VS_f(i),VS_m(j));
            ID(i,j) =  (VS_f(i)>VS_m(j)+eps)-1*(VS_m(j)>VS_f(i)+eps);
            
        end
    end
    
else
    
    error('Value vectors must be columns')
    
end
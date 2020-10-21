%{
Given matrices for 
%}

function U = u_singles(param,C,L,G)

szC = size(C);
szL = size(L);
szG = size(G);

if isequal(szC,szL) && isequal(szL,szG)
    
    if all(C>0) && all(L>0) && all(G>0)
        
        if param.sigma == 1
    
        U = param.sigma_c.*log(C)+param.sigma_l.*log(L)+...
            param.sigma_g.*log(G);
    
        else
    
        U = (param.sigma_c/(1-param.sigma))*C.^(1-param.sigma)+...
        (param.sigma_l/(1-param.sigma))*L.^(1-param.sigma)+...
        (param.sigma_g/(1-param.sigma))*G.^(1-param.sigma);
    
        end
        
    else
       
        message_1 = ['Consumption, leisure and household production '...
            'matrices must not contain non-positive elements'];
        error(message_1);
        
    end    
   
else
    
    message_2 = ['Consumption, leisure and household production '...
        'matrices must have the same dimensions'];
    error(message_2);
    
end


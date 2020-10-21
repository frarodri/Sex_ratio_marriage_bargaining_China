%{
Computes the effective home production time for married households. 
%}

function H = h_married(param,Hf,Hm)

if param.eta == 1
    
    H = (Hf.^param.eta_f).*(Hm.^(1-param.eta_f));
    
else
    
    H = (param.eta_f*Hf.^(1-param.eta)+...
        (1-param.eta_f)*Hm.^(1-param.eta)).^(1/(1-param.eta));
   
end


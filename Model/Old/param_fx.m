%{
Calculates the marriage equilibrium for different values of a parameter
specified by the user
%}

function [Q] = ...
    param_fx(parameter,min,max,step,param,theta_S,OMEGA_f,OMEGA_m,Z_f,Z_m)

values = (min:step:max);

Q = zeros(5,5,length(values));

for i = 1:length(values)
   
     param.parameter = values(i);
    
    [Q(:,:,i), ~, ~, ~, ~] = ...
        marriage_eq(param,theta_S,OMEGA_f,OMEGA_m,Z_f,Z_m);
    
end

end
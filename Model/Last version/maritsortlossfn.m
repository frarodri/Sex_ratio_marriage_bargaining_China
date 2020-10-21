%{
Computes the loss function that will be minimized to select the means for
the match quality draws between agents with different types  
%}

function L = maritsortlossfn(param,targ,theta_S,Q,MU_vec)

types = size(param.Pf,1);
param.MU = vec2mat(MU_vec,types);

MS = marital_sorting(param,theta_S,Q);

L = sum(sum((MS-targ.MS).^2));


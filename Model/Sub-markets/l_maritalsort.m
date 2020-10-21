%{
Computes the loss function that will be minimized to select the means for
the match quality draws between agents with different types  
%}

function L = l_maritalsort(param,targ,theta_S,Q,MU_vec)

types = size(Q,1);

MU_aux = vec2mat(MU_vec,types);
MU_l = [zeros(1,types);tril(MU_aux)];
MU_u = [triu(MU_aux,1);zeros(1,types)];
MU = MU_l+MU_u;
param.MU = MU;

MS = marital_sorting(param,theta_S,Q);

L = sum(sum((MS-targ.MS).^2));
%drawmean_finder

types = size(Q,1);
x0 = zeros(types*(types-1),1);

mus = fminunc(@(MU_vec)l_maritalsort(param,target,theta_S,Q,MU_vec),x0);

MU_aux = vec2mat(mus,types);
MU_l = [zeros(1,types);tril(MU_aux)];
MU_u = [triu(MU_aux,1);zeros(1,types)];
MU = MU_l+MU_u;


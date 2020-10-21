function [lf,lm,hf,hm,nf,nm,mus] = ...
    calib_iteration(param,wages,targ,USf,USm)

types = size(param.Pf,1);

t=cputime;

%% Compute the SS
[theta_S,Q,~] = ...
    sseq(param,wages,USf,USm,'bargaining','Egalitarian');

%% Minimize marital sorting error
x0 = reshape(param.MU,[1,numel(param.MU)]);

options = optimoptions('fminunc','Display','off');
mus = fminunc(@(MU_vec)maritsortlossfn(param,targ,theta_S,Q,MU_vec),...
    x0,options);

%% Recompute the steady state equilibrium with the new means

param.MU = vec2mat(mus,types);
[theta_S,Q,PWf] = sseq(param,wages,USf,USm,'bargaining','Egalitarian');

%% Compute model moments

MS = marital_sorting(param,theta_S,Q);

[~,~,Lf,Lm,~,Hf,Hm,~,~,Nf,Nm] = ...
    per_sol_married(param,wages.f,wages.m,PWf,'quietly','true');
                 
lf = sum(sum(MS.*Lf));
lm = sum(sum(MS.*Lm));
hf = sum(sum(MS.*Hf));
hm = sum(sum(MS.*Hm));
nf = sum(sum(MS.*Nf));
nm = sum(sum(MS.*Nm));

%% Check how far from targets we are

% ta_results.Statistic = {'Wives housework';'Wives paid work';...
%     'Wives leisure';'Husbands housework';'Husbands paid work';...
%     'Husbands leisure'};
% 
% ta_results.Data = [targ.hf;targ.nf;targ.lf;targ.hm;targ.nm;targ.lm];
% 
% ta_results.Model = [hf;nf;lf;hm;nm;lm];
% 
% Time_allocation_results = struct2table(ta_results)

%% Compute the loss

% L = (lf-targ.lf)^2+(lm-targ.lm)^2+(hf-targ.hf)^2+(hm-targ.hm)^2+...
%     (nf-targ.nf)^2+(nm-targ.nm)^2+sum(sum((MS-targ.MS).^2));

e = cputime-t;

message_elapsed_time = ['Elapsed time for one calibration iteration: %g seconds \n'];
fprintf(message_elapsed_time,e);
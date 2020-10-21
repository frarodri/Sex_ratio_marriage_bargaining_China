%% Compute the PREE with Nash Bargaining in a submarket

pwf = 1/2;

qr = fzero(@(qr)mgains_initial(param,wages,theta_S,pwf,mu,qr),0);

[USf,USm] = uflow_singles(param,wages);

[UMf,UMm] = uflow_married(param,wages.f,wages.m,pwf);

[VSf, VSm] = ...
    v_singles_submkts_initial(param,theta_S,UMf,UMm,USf,USm,qr,mu);

% Set initial value for error and tolerance for loop
current_error = 1;
tol = 10^(-9);
it = 1;

%% Loop 
t = cputime;
while current_error>tol
   
   pwf
   qr
   qr_prior = qr;
   pwf_prior = pwf;
   
   start_point = pwf;
    
   [EVMf,EVMm,pwf] = VM_expect(param,wages,VSf,VSm,qr,mu,start_point);
   [VSf, VSm] = v_singles_submkts(param,theta_S,USf,USm,EVMf,EVMm,qr,mu);
   
   qr = fzero(@(q)mgains(param,wages,VSf,VSm,pwf,q),0);
   
   error_report = max(abs(qr-qr_prior),abs(pwf-pwf_prior));
    
%     if strcmp(p.Results.quietly,"false")
        msg = ['
         fprintf('Iteration %d, Pareto current error is %g\n',it,error_report);

%     end

    current_error = error_report;
    it = it+1;
    
end
e = cputime-t;
fprintf('Loop elapsed time: %g\n seconds',e);
%%
% Check

Gpwf = linspace(0.01,0.99);
Q = zeros(size(Gpwf));

for i=1:length(Gpwf)
    
    pwf = Gpwf(i);
    
    q_r = fzero(@(q)mgains(param,wages,VSf,VSm,pwf,q),0);
    
    Q(i) = q_r;
    
end
    
plot(Gpwf,Q);

%%


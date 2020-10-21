function [Gf,Gm,pwf] = ...
    nashbarg(param,wage_f,wage_m,VSf,VSm,qr,mu,start_point)

%% Attempt to do it using Matlab's provided integration function
% EVMf = integral(@(q)integrand_nash_f(param,wages,VSf,VSm,mu,q),qr,inf);
% EVMm = integral(@(q)integrand_nash_m(param,wages,VSf,VSm,mu,q),qr,inf);

%% Doing it brute force
npoints = 100;
Grid = linspace(qr,qr+3,npoints);

binsz = 3/(npoints-1);

Gmidpts = (Grid(:,1:end-1)+Grid(:,2:end))/2;

Intf = zeros(size(Gmidpts));
Intm = zeros(size(Gmidpts));

PWf = zeros(size(Gmidpts));

for i=1:length(Gmidpts)
    
    q = Gmidpts(i);
    pwf = nashbarg_pw(param,wage_f,wage_m,VSf,VSm,q,...
        'start_point',start_point);
    
    [UMf,UMm] = uflow_married(param,wage_f,wage_m,pwf,'quietly','true');
    Mf = (UMf+q)/(1-param.beta*(1-param.delta));
    Mm = (UMm+q)/(1-param.beta*(1-param.delta));
    
    PWf(i) = pwf;
    
    Intf(i) = Mf*normpdf(q,mu);
    Intm(i) = Mm*normpdf(q,mu);
    
end

Gf = sum(Intf)*binsz;
Gm = sum(Intm)*binsz;

pwf = PWf(1);
%uflow_married(param,wage_f,wage_m,pwf);

% plot(Gmidpts,PWf);   

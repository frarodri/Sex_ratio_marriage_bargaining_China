function [EVMf,EVMm,pwf] = VM_expect(param,wages,VSf,VSm,qr,mu,start_point)

%% Attempt to do it using Matlab's provided integration function
% EVMf = integral(@(q)integrand_nash_f(param,wages,VSf,VSm,mu,q),qr,inf);
% EVMm = integral(@(q)integrand_nash_m(param,wages,VSf,VSm,mu,q),qr,inf);

%% Doing it brute force
npoints = 100;
G = linspace(qr,qr+3,npoints);

binsz = 3/(npoints-1);

Gmidpts = (G(:,1:end-1)+G(:,2:end))/2;

Intf = zeros(size(Gmidpts));
Intm = zeros(size(Gmidpts));

PWf = zeros(size(Gmidpts));

for i=1:length(Gmidpts)
    
%     q = Gmidpts(i);
%     Intf(i) = integrand_nash_f(param,wages,VSf,VSm,mu,q);
%     Intm(i) = integrand_nash_m(param,wages,VSf,VSm,mu,q);
%     
%     PWf(i) = nashbarg_pw(param,wages,VSf,VSm,q);

    q = Gmidpts(i);
    pwf = nashbarg_pw(param,wages,VSf,VSm,q,'start_point',start_point);
    
    [UMf,UMm] = uflow_married(param,wages.f,wages.m,pwf,'quietly','true');
    Mf = (UMf+q)/(1-param.beta*(1-param.delta));
    Mm = (UMm+q)/(1-param.beta*(1-param.delta));
    
    PWf(i) = pwf;
    
    Intf(i) = Mf*normpdf(q,mu);
    Intm(i) = Mm*normpdf(q,mu);
    
end

EVMf = sum(Intf)*binsz;
EVMm = sum(Intm)*binsz;

pwf = PWf(1);

plot(Gmidpts,PWf);   

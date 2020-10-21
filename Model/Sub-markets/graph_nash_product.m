%%
% Graph the Nash product as a function of the Pareto weight of the wife and
% the sex ratio

Thetas = [0.5 ; 1 ; 2 ; 3];
szThetas = length(Thetas);

Chis = linspace(0.01,0.99);
szChis = length(Chis);

NP = zeros(szThetas,szChis);

for i = 1:length(Thetas)

    theta_S = Thetas(i);

    for j = 1:length(Chis)

    pwf = Chis(j);
    np = nash_product(param,wages,theta_S,qr,pwf,mu);
    NP(i,j) = np;
    
    end
    
end

legend_text = "$$\theta_S$$"+"="+Thetas;
title_x_axis = "$$\chi_f$$";

plot(Chis,NP);
legend(legend_text,'Interpreter','latex','Location','best')
legend('boxoff')
xlabel(title_x_axis,'interpreter','latex')
ylabel('Nash product')
    
%%
% Graph the utility flow of married couples as a function of the Pareto 
% weight

Chis = linspace(0.01,0.99);
szChis = length(Chis);

UMf = zeros(size(Chis));
UMm = zeros(size(Chis));

for i=1:szChis
    
    pwf = Chis(i);
    
    [um_f,um_m] = ...
        uflow_married(param,wages.f,wages.m,pwf,'quietly','true');
    
    UMf(i) = um_f;
    UMm(i) = um_m;
    
end

plot(Chis,UMf,Chis,UMm);

%% 
% Graph the gains from marriage as a function of the Pareto weight of the
% wife

Chis = linspace(0.01,0.99);
szChis = length(Chis);

Wf = zeros(size(Chis));
Wm = zeros(size(Chis));

for i=1:szChis
    
    pwf = Chis(i);
    
    % Compute the flow value of being single
    [USf,USm] = uflow_singles(param,wages);

    % Compute the flow value of marriage with current Pareto weight for wife
    [UMf,UMm] = uflow_married(param,wages.f,wages.m,pwf,'quietly','true');

    % Compute the value of being single with current Pareto weight for wife
    [VSf, VSm] = v_singles_submkts(param,theta_S,UMf,UMm,USf,USm,qr,mu);

    % Compute the gains from marriage for wife and husband 
    Wf(i) = (UMf+qr)/(1-param.beta*(1-param.delta))-VSf;
    Wm(i) = (UMm+qr)/(1-param.beta*(1-param.delta))-VSm;
    
end

plot(Chis,Wf,Chis,Wm);

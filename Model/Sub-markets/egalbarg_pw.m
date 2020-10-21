%{
Computes the Pareto weight for the wife that implements the Egalitarian 
Bargaining solution within marriage, i.e. the Pareto weight that equalizes
the gains from marriage for both spouses, given a reservation quality
for marriage.
%}

function pw = egalbarg_pw(param,theta_S,wage_f,wage_m,qr,mu)

% Set wages and expectations
wages.f = wage_f;
wages.m = wage_m;

% A = [];
% b = []; 
% Aeq = [];
% beq = [];
% x0 = rand();

% Old: uses fmincon, caused problems
%options = optimoptions('fmincon','Display','off');
%pw = fmincon(@(pwf)mgains_diff(param,wages,theta_S,qr,pwf,mu)^2,...
%                                x0,A,b,Aeq,beq,0,1,[],options);

% New: uses fzero
pw = ...
    fzero(@(pwf)mgains_diff(param,wages,theta_S,qr,pwf,mu),[10^-(10),1]);



        
        
        
        
        
        





    
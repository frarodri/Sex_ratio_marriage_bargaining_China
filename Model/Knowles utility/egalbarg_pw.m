%{
Computes the Pareto weight for the wife that implements the Egalitarian 
Bargaining solution within marriage, i.e. the Pareto weight that equalizes
the gains from marriage for both spouses, given a reservation quality
for marriage. It uses fmincon to minimize the difference in the gains of 
marriage.
%}

function pw = egalbarg_pw(param,theta,wage_f,wage_m,q)

% Set wages and expectations
wages.f = wage_f;
wages.m = wage_m;

expect.Pf = 1;
expect.Pm = 1;
expect.theta_S = theta;

A = [];
b = []; 
Aeq = [];
beq = [];
x0 = rand();

options = optimoptions('fmincon','Display','off');
pw = fmincon(@(mu)mgains_diff(param,expect,wages,q,mu)^2,...
                                x0,A,b,Aeq,beq,0,1,[],options);
        
        
        
        
        
        





    
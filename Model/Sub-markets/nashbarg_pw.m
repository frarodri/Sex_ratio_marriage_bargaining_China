%{
Computes the Pareto weight for the wife that implements the Nash 
Bargaining solution within marriage, i.e. the Pareto weight that maximizes 
the Nash product , given a reservation quality for marriage. 
%}

function pw = nashbarg_pw(param,theta_S,wage_f,wage_m,qr,mu)

% Set wages 
wages.f = wage_f;
wages.m = wage_m;

A = [];
b = []; 
Aeq = [];
beq = [];
x0 = rand();

options = optimoptions('fmincon','Display','off');
pw = fmincon(@(pwf)-nash_product(param,wages,theta_S,qr,pwf,mu),...
                               x0,A,b,Aeq,beq,0,1,[],options);
                           
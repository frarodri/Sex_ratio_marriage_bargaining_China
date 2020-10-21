%{
Computes the Pareto weight for the wife that implements the Nash 
Bargaining solution within marriage, i.e. the Pareto weight that maximizes 
the Nash product , given a reservation quality for marriage. 
%}

function pw = nashbarg_pw_entrants(param,vs_f,vs_m,wage_f,wage_m,qr)


% Set wages 
wages.f = wage_f;
wages.m = wage_m;

A = [];
b = []; 
Aeq = [];
beq = [];
x0 = rand();

options = optimoptions('fmincon','Display','off');
pw = fmincon(@(pwf)-nash_product_entrants(param,wages,...
    vs_f,vs_m,qr,pwf),x0,A,b,Aeq,beq,0,1,[],options);
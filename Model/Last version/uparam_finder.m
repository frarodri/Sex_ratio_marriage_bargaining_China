

x0 = 1/3*[1 1 1];
A = [];
b = []; 
Aeq = [1 1 1];
beq = 1;
lb = [0 0 0];
ub = [];


options = optimoptions('fmincon','Display','off');
sigmas = fmincon(@(x)lossfn_married(param,target,wage_f,wage_m,pwf,x),...
                                        x0,A,b,Aeq,beq,lb,ub,[],options);
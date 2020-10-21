

param.phi = 0.75;
param.w_m = 1;
param.w_f = param.phi*param.w_m;
param.gamma = 1.9;

data.growth = 7;
data.l_init = 0.89;
data.l_final = 0.74;

A = [];
b = []; 
Aeq = [];
beq = [];

x0 = [0.5,0.5];

lb = [0,0];

ub = [1,0.7];


calibrated_params = ...
    fmincon(@(x)cost_function_test(data,param,x),x0,A,b,Aeq,beq,lb,ub);

lambdas = (0:0.1:1);
cbars = (0:0.1:0.7);

cost = zeros(length(lambdas),length(cbars));

for i = 1:length(lambdas)
    for j = 1:length(cbars)
        
        x = [lambdas(i),cbars(j)];
        cost(i,j) = cost_function_test(data,param,x);
        
    end
end

surf(cbars,lambdas,cost)
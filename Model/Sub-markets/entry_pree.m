%{
Computes the Partial Rational Expectations Equilibrium (PREE) with
Egalitarian Bargaining for the entrants marriage market. That is, given 
continuation values for all types of agents, i.e. the value of being 
single in each of the same-type submarkets, computes the reservation match 
qualities and Pareto weights for marriage between entrants that are 
presented with a marriage opportunity. Stores the results in matrices Q and 
PW, respectively.  
%}

function [Q,PWf] = entry_pree(param,wages,VSf,VSm,varargin)

p = inputParser;
paramName = 'quietly';
defaultVal = "false";
errorMsg = 'Quietly must be either true or false';
validationFcn = @(x) assert(strcmp(x,"true") | strcmp(x,"false"),errorMsg);
addParameter(p,paramName,defaultVal,validationFcn);
parse(p,varargin{:});

szf = size(param.Pf,1);
szm = size(param.Pm,1);

VSf = VSf.*ones(szf,szm);
VSm = VSm'.*ones(szf,szm);

%%
% % Egalitarian Bargaining Pareto weights
% 
% % First I need to define the Pareto weights. I will use the bisection
% % method, because the difference in gains from marriage is strictly
% % increasing in this case since the value of being single does not depend
% % on the Pareto weights in this stage. Also, it is easier to implement with
% % matrices than fzero
% 
% % Start by defining the initial Pareto weights and upper and lower bounds
% PWf = 0.5*ones(szf,szm);
% UB = ones(szf,szm);
% LB = zeros(szf,szm);
% 
% % Compute the flow value of marriage with current Pareto weight for wife
% [UMf,UMm] = uflow_married(param,wages.f,wages.m,PWf,'quietly','true');
% 
% % Compute the gains from marriage and the difference with current Pareto 
% % weight for wife
% Wf = UMf/(1-param.beta*(1-param.delta))-VSf;
% Wm = UMm/(1-param.beta*(1-param.delta))-VSm;
% DW = Wf-Wm;
% 
% % Set tolerance level for convergence
% tol = 10^(-5);
% 
% % Start bisection algortihm
% 
% while sum(sum(abs(DW)))>tol
%     
%     UPD=(DW>0);
%     
%     % Update Pareto weights and upper and lower bounds
%     UB = UPD.*PWf+(1-UPD).*UB;
%     LB = (1-UPD).*PWf+UPD.*LB;
%     PWf = (UB+LB)/2;
%     
%     % Re-compute the flow value of marriage with current Pareto weight
%     % for wife
%     [UMf,UMm] = uflow_married(param,wages.f,wages.m,PWf,'quietly','true');
%     
%     % Re-compute the gains from marriage and the difference with 
%     % current Pareto weight for wife
%     Wf = UMf/(1-param.beta*(1-param.delta))-VSf;
%     Wm = UMm/(1-param.beta*(1-param.delta))-VSm;
%     DW = Wf-Wm;
%     
% end
% 
% % Compute the reservation quality of marriage by setting it to the value
% % that equalizes the value of marriage and the value of being single
% Qf = (1-param.beta*(1-param.delta))*VSf-UMf;
% Qm = (1-param.beta*(1-param.delta))*VSm-UMm;
% Q = (Qf+Qm)/2;

%%
% Nash bargaining Pareto weights

Q = zeros(szf,szm);

for f=1:szf
    
    wage_f = wages.f(f);
    
    for m=1:szm
        
        wage_m = wages.m(m);
        
        vs_f = VSf(f,m);
        vs_m = VSm(f,m);
        
        % Guess an initial reservation match qualitiy q
        qr = 0;
        
        % Set initial value for error and tolerance for loop
        current_error = 1;
        tol = 10^(-5);
        it = 1;
        
        while current_error>tol
   
            qr_prior = qr;
        
            %% 
            % 1. Compute the Pareto weights that implement the Nash 
            % Bargaining solution for the current reservation match quality

            pw = nashbarg_pw_entrants(param,vs_f,vs_m,wage_f,wage_m,qr);
            
            PWf(f,m) = pw;
            
            %% 
            % 2. Compute the flow value of marriage for the agents with the
            % current Pareto weights
            
            [um_f,um_m] = uflow_married(param,wage_f,wage_m,pw,...
                'quietly','true');
            
            %% 
            % 3. Update the reservation match quality to set the gains from
            % marriage to zero
            
            qr = ...
                ((1-param.beta*(1-param.delta))*(vs_f+vs_m)-(um_f+um_m))/2;
            
            %%
            % 4. Update the error
            
            current_error = abs(qr-qr_prior);
            
            
        end
        
        Q(f,m) = qr;
        
    end
    
end

    
    
    
    
    
    

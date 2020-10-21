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
paramName1 = 'quietly';
defaultVal1 = "false";
errorMsg1 = 'Quietly must be either true or false';
validationFcn1 = ...
    @(x) assert(strcmp(x,"true") | strcmp(x,"false"),errorMsg1);
addParameter(p,paramName1,defaultVal1,validationFcn1);

paramName2 = 'bargaining';
defaultVal2 = "Nash";
errorMsg2 = 'The bargaining solution must be Egalitarian or Nash';
validationFcn2 = ...
    @(x) assert(strcmp(x,"Egalitarian") | strcmp(x,"Nash"),errorMsg2);
addParameter(p,paramName2,defaultVal2,validationFcn2);

parse(p,varargin{:});

szf = size(param.Pf,1);
szm = size(param.Pm,1);

VSf = VSf.*ones(szf,szm);
VSm = VSm'.*ones(szf,szm);

if strcmp(p.Results.bargaining,"Nash")
    
    % Nash bargaining
    
    Q = zeros(szf,szm); 
    PWf = zeros(szf,szm);

    for f=1:szf

        wage_f = wages.f(f);

        for m=1:szm

            wage_m = wages.m(m);

            vs_f = VSf(f,m);
            vs_m = VSm(f,m);

            % Set initial Pareto weight of wife 
            pwf = 1/2;

            % Compute the initial reservation match quality 
            qr = fzero(@(q)mgains(param,wage_f,wage_m,vs_f,vs_m,pwf,q),0);

            % Set initial value for error and tolerance for loop
            current_error = 1;
            tol = 10^(-5);
            it = 1;

            while current_error>tol

                qr_prior = qr;
                pwf_prior = pwf;
                
                % 1. Compute the Pareto weights that implement the Nash 
                % Bargaining solution for the current reservation match 
                % quality
                pwf = nashbarg_pw(param,wage_f,wage_m,vs_f,vs_m,qr,...
                    'start_point',pwf_prior);

                PWf(f,m) = pwf;

                %% 
                % 2. Update the reservation match quality to set the gains 
                % from marriage to zero for the current Pareto weights

                qr = fzero(...
                    @(q)mgains(param,wage_f,wage_m,vs_f,vs_m,pwf,q),...
                    qr_prior);
                Q(f,m) = qr;

                %%
                % 3. Update the error

                error_report = max(abs(qr-qr_prior),abs(pwf-pwf_prior));

                if strcmp(p.Results.quietly,"false")

                    message = 'Iteration %d, current error is %g\n';

                     fprintf(message,it,error_report);

                end

                current_error = error_report;
                it = it+1;

            end

        end

    end
    
else
    
    % Egalitarian Bargaining 
    
    %{
    I use the bisection method, since the difference in gains from mariage
    is strictly increasing for entrants since the value of being single is
    known, and the implementation is parsimonious than using fzero with a
    loop that runs through every potential type-pair
    %}

    % Set the initial Pareto weights and upper and lower bounds
    PWf = (1/2)*ones(szf,szm);
    ub = 1-10^(-10);
    UB = ones(szf,szm)*ub;
    lb = 10^(-10);
    LB = ones(szf,szm)*lb;

    % Compute the flow value of marriage with current Pareto weight for wife
    [UMf,UMm] = uflow_married(param,wages.f,wages.m,PWf,'quietly','true');

    % Compute the gains from marriage and the difference with current Pareto 
    % weight for wife
    Wf = UMf/(1-param.beta*(1-param.delta))-VSf;
    Wm = UMm/(1-param.beta*(1-param.delta))-VSm;
    DW = Wf-Wm;

    % Set tolerance level for convergence
    tol = 10^(-5);
    PWf_diff = 1;

    % Start bisection algortihm

    while sum(abs(DW),'all')>tol && sum(abs(PWf_diff),'all')>tol

        UPD=(DW>0);
        PWf_prior = PWf;
        
        % Update Pareto weights and upper and lower bounds
        UB = UPD.*PWf+(1-UPD).*UB;
        LB = (1-UPD).*PWf+UPD.*LB;
        PWf = (UB+LB)/2;

        % Re-compute the flow value of marriage with current Pareto weight
        % for wife
        [UMf,UMm] = ...
            uflow_married(param,wages.f,wages.m,PWf,'quietly','true');

        % Re-compute the gains from marriage and the difference with 
        % current Pareto weight for wife
        Wf = UMf/(1-param.beta*(1-param.delta))-VSf;
        Wm = UMm/(1-param.beta*(1-param.delta))-VSm;
        DW = Wf-Wm;
        PWf_diff = PWf-PWf_prior;

    end

    % Compute the reservation quality of marriage by setting it to the value
    % that equalizes the value of marriage and the value of being single
    Qf = (1-param.beta*(1-param.delta))*VSf-UMf;
    Qm = (1-param.beta*(1-param.delta))*VSm-UMm;
    Q = (Qf+Qm)/2;
    
end


%%


%%



    
    
    
    
    
    

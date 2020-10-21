%{
Given vectors of different wages for females and males and Pareto 
weights for females, computes whenever possible the analytical solution 
to the problem faced each period by married agents. If the analytical 
solution is not valid, computes a numerical one.
%}

function [Cf,Cm,Lf,Lm,G,Hf,Hm,H,EQ,Nf,Nm] = ...
                      per_sol_married(param,wages_f,wages_m,PWf,varargin)

p = inputParser;
paramName = 'quietly';
defaultVal = 'false';
errorMsg = 'Quietly must be either true or false';
validationFcn = @(x) assert(strcmp(x,"true") | strcmp(x,"false"),errorMsg);
addParameter(p,paramName,defaultVal,validationFcn);
parse(p,varargin{:});
                              
if (iscolumn(wages_f) && iscolumn(wages_m)) && ...
        (all(wages_f>0) && all(wages_m>0))
    
    szf = size(wages_f);
    szm = size(wages_m);
           
    if isequal(size(PWf),[szf(1),szm(1)]) && all(all(PWf>0)) && ...
            all(all(PWf<1))
        
        Wf = wages_f.*ones(szf(1),szm(1));
        Wm = wages_m'.*ones(szf(1),szm(1));
        PWm = 1-PWf;

        XF = (param.eta_f*Wm./((1-param.eta_f)*Wf)).^(1/param.eta);
        XH = (param.eta_f.*XF.^(1-param.eta)+1-param.eta_f).^...
             (1/(1-param.eta));
        XE = ((1-param.alpha_g).*Wm.*XH.^(1-param.eta))./...
             (param.alpha_g*param.pe*(1-param.eta_f));
        XG = param.A_g*(XE.^(1-param.alpha_g)).*(XH.^param.alpha_g);
        D = (param.pe*(XE./XH).^param.alpha_g)./...
            (param.A_g*(1-param.alpha_g));

        NUM_LM = (PWf*param.sigma_c).^(1/param.sigma)+...
                 (PWm*param.sigma_c).^(1/param.sigma)+...
                 Wf.*(PWf*param.sigma_l./Wf).^(1/param.sigma)+...
                 Wm.*(PWm*param.sigma_l./Wm).^(1/param.sigma)+...
                 ((param.sigma_g./D).^(1/param.sigma)).*...
                 (Wf.*(XF./XG)+Wm./XG+param.pe*(XE./XG));
        LM = (NUM_LM./(Wf+Wm)).^param.sigma;

        Cf = (PWf*param.sigma_c./LM).^(1/param.sigma);
        Cm = (PWm*param.sigma_c./LM).^(1/param.sigma);
        G = (param.sigma_g./(LM.*D)).^(1/param.sigma);
        Hf = G.*XF./XG;
        Hm = G./XG;
        Lf = (PWf*param.sigma_l./(LM.*Wf)).^(1/param.sigma);
        Lm = (PWm*param.sigma_l./(LM.*Wm)).^(1/param.sigma); 
        H = h_married(param,Hf,Hm);
        EQ = XE.*G./XG;
        Nf = 1-Lf-Hf;
        Nm = 1-Lm-Hm;

        % Analytical solution only works when time constraints are not 
        % binding, we need to check whether Nf, Nm, Hf or Hm contain 
        % negative elements and substitute the solution for a numerical 
        % one if that is the case, or if the solution for any of the above 
        % quantities is not a number (NaN)
        
        NaN_check = Cf+Cm+G+Hf+Hm+Lf+Lm+H+EQ+Nf+Nm;

        if all(all(Nf>=0)) && all(all(Nm>=0)) && all(all(Hf>=0)) && ...
           all(all(Hm>=0)) && all(all(~isnan(NaN_check)))
       
            if strcmp(p.Results.quietly,"false")
                
                disp('All solutions are interior');
                
            end

        else
            
            if strcmp(p.Results.quietly,"false")
                
                disp('Some solutions not interior');
                
            end

            CS = (Nf<0)+(Nm<0)+(Hf<0)+(Hm<0)+isnan(NaN_check);

            for i=1:size(CS,1)

                for j=1:size(CS,2)

                    if CS(i,j)~=0

                        wf = Wf(i,j);
                        wm = Wm(i,j);
                        pwf = PWf(i,j);

                        A = [0,0,1,0,1,0,0; 0,0,0,1,0,1,0;...
                             1,1,wf,wm,wf,wm,param.pe];
                        b = [1;1;wf+wm];
                        Aeq = [];
                        beq = [];
                        x0 = [0.2,0.2,0.1,0.1,0.1,0.1,0.1];

                        lb = [0,0,0,0,0,0,0];
                        ub = [Inf,Inf,1,1,1,1,Inf];
                        options = optimoptions('fmincon','Display','off');

                        sol = ...
                            fmincon(@(ctrl)-obj_married(param,ctrl,pwf),...
                                      x0,A,b,Aeq,beq,lb,ub,[],options);

                        cf = sol(1);
                        cm = sol(2);
                        hf = sol(3);
                        hm = sol(4);                
                        lf = sol(5);
                        lm = sol(6);
                        h = h_married(param,hf,hm);
                        eq = sol(7);
                        g = home_prod(param,h,eq);
                        nf = 1-hf-lf;
                        nm = 1-hm-lm;

                        Cf(i,j) = cf;
                        Cm(i,j) = cm;
                        G(i,j) = g;
                        Hf(i,j) = hf;
                        Hm(i,j) = hm;
                        Lf(i,j) = lf;
                        Lm(i,j) = lm; 
                        H(i,j) = h;
                        EQ(i,j) = eq;
                        Nf(i,j) = nf;
                        Nm(i,j) = nm;

                    end

                end

            end

        end
        
    else
        
        message1 = ['Matrix of Pareto weights for wife must have equal '...
            'number of rows as wage vector for females, equal '...
            'number of columns as wage vector for males and all of its '...
            'entries should be between 0 and 1'];
        error(message1);
        
    end
    
else
    
    message2 = ['Wage vectors must be columns without non-positive '...
                'elements'];
    error(message2);
    
end
                                    
                                    

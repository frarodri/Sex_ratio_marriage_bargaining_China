function [Q,theta_S,P_f,P_m] = marriage_free(param,wages,init_expect)

szP_f = size(param.P_f,1);
szP_m = size(param.P_m,1);

szOMEGA_f = size(wages.OMEGA_f,1);
szOMEGA_m = size(wages.OMEGA_m,1);

if iscolumn(param.P_f) && iscolumn(param.P_m) && ...
        iscolumn(wages.OMEGA_m) && iscolumn(wages.OMEGA_m)
    
    if szP_f==szOMEGA_f && szP_m==szOMEGA_m
        
        expect.theta_S = init_expect.theta_S;
        expect.P_f = init_expect.P_f;
        expect.P_m = init_expect.P_m;
        
        err = 1; % Initialize error
        tol = 1*exp(-10); % Tolerance
        
        iter = 1;    

        while err>tol
               
            [Q,~,~,~,~] = marriage_pree(param,expect,wages);
            
            [theta_S, P_f, P_m] = singles_dist(Q,param);
                       
            err = abs(theta_S-expect.theta_S)+sum(abs(P_f-expect.P_f))+...
                  sum(abs(P_m-expect.P_m)); 
            
            expect.theta_S = theta_S;
            expect.P_f = P_f;
            expect.P_m = P_m;  
              
            iter = iter+1;
            
        end
        
    else
    
        error_text_1 = ['Distribution of types among new entrants and '...
            'wage vectors must have the same dimensions for each sex'];
        error(error_text_1);
        
    end
    
else
    
    error_text_2 = ['Distribution of types among new entrants and '...
        'wages must be vectors'];
    error(error_text_2);
    
end

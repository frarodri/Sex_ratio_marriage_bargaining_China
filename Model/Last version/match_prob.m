%{
Calculates the matching probabilities for females and males given 
matching efficiency and elasticity parameters and the sex ratio among 
singles 
%}

function [PIf, PIm] = match_prob(param,theta_S)

p = inputParser;
errormsgtheta_S = ['The sex ratio must be a scalar or vector of strictly '...
    'positive numbers.'];
valfntheta_S = @(x) assert(isvector(x) && min(x>0),errormsgtheta_S);
addRequired(p,'theta_S',valfntheta_S);
parse(p,theta_S);

if iscolumn(theta_S)
    
    Xf = [param.A_x*theta_S.^param.alpha_x theta_S ones(size(theta_S))];

    Xm = [param.A_x*(1./theta_S).^(1-param.alpha_x) ones(size(theta_S)) ...
    1./theta_S];
    
    PIf = min(Xf,[],2);
    PIm = min(Xm,[],2);
    
else
    
    Xf = [param.A_x*theta_S.^param.alpha_x;theta_S;ones(size(theta_S))];

    Xm = [param.A_x*(1./theta_S).^(1-param.alpha_x);ones(size(theta_S));...
    1./theta_S];
    
    PIf = min(Xf);
    PIm = min(Xm);
    
end

% if isscalar(param.A_x) && param.A_x>=0
%     
%     if isscalar(param.alpha_x) && param.alpha_x>=0 && param.alpha_x<=1
%         
%         if isscalar(theta_S) && theta_S>0
%             
%             X_m = [param.A_x*(1/theta_S)^(1-param.alpha_x) 1 1/theta_S];
%             X_f = [param.A_x*theta_S^param.alpha_x theta_S 1];
%             PIf = min(X_f);
%             PIm = min(X_m);
%             
%         else
%             
%             error_text_1 = ['Male to female ratio among singles must '...
%                 'be a positive scalar'];
%             error(error_text_1)
%             
%         end
%         
%     else
%         
%         error_text_2 = ['Matching function elasticity must be a scalar '...
%             'between zero and one'];
%         error(error_text_2)
%         
%     end
%     
% else
%     
%     error_text_3 = ['Matching function efficiency parameter must be a '...
%         'non-negative scalar'];
%     error(error_text_3)
%     
% end


   
    


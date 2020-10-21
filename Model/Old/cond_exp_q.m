%{
Computes the expected quality of matches that turn into marriage given a 
matrix of quality cut-off values Q and records it in a matrix where rows 
are wife's type and columns husband's type, when the distribution of 
the quality of match draws is normal
%}

function EqM = cond_exp_q(param,Q)

szQ = size(Q);
szMU = size(param.MU);

if isequal(szQ,szMU) 
    
    EqM = param.MU+normpdf(Q-param.MU)./(1-normcdf(Q-param.MU));
    
    %{
    Due to Mat-lab's numerical approximation error when calculating
    the standard normal probability and cumulative density 
    functions, the formula above yields either values that are below the 
    cutoff or infinity for large cut-offs. We have to fix this.
    %}
    
    % Conditional expectation cannot be smaller than cut-off, so set it 
    % to cut-off value whenever that's the case 
    EqM = max(EqM,Q);
    
    % Set infinite values to cut-off as well
    Q_inf = Q.*isinf(EqM);
    EqM(isinf(EqM))=0; 
    EqM = EqM+Q_inf;
   
else
    
    error_text = ['Cut-off quality matrix and matrix of means must '...
        'have the same dimensions'];
    error(error_text)

end

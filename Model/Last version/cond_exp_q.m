%{
Computes the expected quality of matches that turn into marriage given a 
matrix of reservation quality values Q and records it in a matrix where 
rows are wife's type and columns husband's type, when the distribution of 
the quality of match draws is normal
%}

function EqM = cond_exp_q(Q,MU)

szQ = size(Q);
szMU = size(MU);

if isequal(szQ,szMU) 
    
    EqM = MU+normpdf(Q-MU)./(1-normcdf(Q-MU));
    
    %{
    Due to Mat-lab's numerical approximation error when calculating
    the standard normal probability and cumulative density 
    functions, the formula above sometimes yields values that are either 
    below the reservation quality or infinity for large reservation quality 
    values. We have to fix this.
    %}
    
    % Conditional expectation cannot be smaller than reservation quality, 
    % so set it to reservation quality value whenever that's the case 
    EqM = max(EqM,Q);
    
    % Set infinite values to reservation quality values as well
    Q_inf = Q.*isinf(EqM);
    EqM(isinf(EqM))=0; 
    EqM = EqM+Q_inf;
   
else
    
    error_text = ['Reservation quality and means matrices must '...
        'have the same dimensions'];
    error(error_text)

end

function [delta, h_f, h_m] = sort_cont_mat_2(MS,Z_f,Z_m)

%{
For a marital sorting matrix, compute the delta measure of assortative 
matching, i.e. the trace divided by the trace of the matrix with random 
sorting, and measures of female and male hypergamy, i.e. the sum of the 
entries in the upper (lower) triangle (excluding the main diagonal) 
divided by the same number for a matrix with random sorting 
%}

[rows,cols] = size(MS);

if rows==cols && rows>=2
    
    if rows==size(Z_f,1) && cols==size(Z_m,1)
        
        RAND = Z_f*Z_m';
        
        delta = trace(MS)/trace(RAND);
        h_f = sum(sum(triu(MS,1)))/sum(sum(triu(RAND,1)));
        h_m = sum(sum(tril(MS,-1)))/sum(sum(tril(RAND,-1)));
        
    else
        
        error('Marginal distributions must conform with sorting matrix')
        
    end
           
else
    
    error('Marital sorting matrix must be square and at least 2x2')
    
end


   

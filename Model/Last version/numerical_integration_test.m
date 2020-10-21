% Expected value of a normally distributed random variable, computed
% numerically and anallytically

%% Parameters

% Mean and standard deviation
mu = 0;
sigma = 1;

% Truncation points, a<b
a = 0;
b = inf;

% Auxiliary variables
alpha = (a-mu)/sigma;
beta = (b-mu)/sigma;
Z = normcdf(beta)-normcdf(alpha);

%% Analytical solution

mean_analyt = mu+((normpdf(alpha)-normpdf(beta))/Z)*sigma;

%% Numerical solution

% Parameters for numerical approximation
% Number of points in the grid, the more points, the more accurate the
% approximation will be
npoints = 200;
% Number of standard deviations to extend the grid when one of the
% truncation points is infinity or minus infinity
standard_devs = 3;

% Choose the upper and lower values for the grid
if a~=-inf
    
    l=a;

elseif a==-inf && b~=inf
    
    l=b-standard_devs*sigma;
    
else
    
    l=mu-standard_devs*sigma;
    
end

if b~=inf
    
    u=b;
    
elseif b==inf
    
    u=a+3;
    
else
    
    u=mu+standard_devs*sigma;
    
end

binsz = (u-l)/(npoints-1);

% Create a grid of points
G = linspace(l,u,npoints);
Gmidpts = (G(:,1:end-1)+G(:,2:end))/2;

% Auxiliary variables (vectors)
Xi = (Gmidpts-mu)/sigma;

% Compute the value of the pdf for each point in the grid
Gpdf = normpdf(Xi)/(Z*sigma);

% Compute the numerical approximation
mean_num = sum((Gmidpts.*Gpdf)*binsz);

% Approximation error
error = abs(mean_analyt-mean_num); 
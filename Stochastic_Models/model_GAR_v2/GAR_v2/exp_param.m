%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This function computes the parameter "a" of an exponential fitting:
%
%                           yfit = exp(-x/a)   
%
% Input:
% x : vector of the independent observed data.
% y : vector of the dependent observed data.
%
% Output:
% a : parameter of the exponential model.
%
% Christian González - PUC. June/2016.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [a] = exp_param(x,y)

% Convert to column vector
x = x(:);
y = y(:);

% Find a convenient starting point  
i = find(y<=0); i=min(i);
if isempty(i)
    mu = sum(x.*y)/sum(y);
    A = -x./(log(trapz(x,y)/mu)-x/mu);
    sp = A(end);    
else
    sp = x(i);
	%i=find(y==min(y)); i=min(i);
end

% Compute the parameter "a" of the exponential curve with a Least-Squares Fitting
f                  = fittype('exp(-x/a)');
[fit1,gof,fitinfo] = fit(x,y,f,'StartPoint',sp);
a                  = fit1.a;

return
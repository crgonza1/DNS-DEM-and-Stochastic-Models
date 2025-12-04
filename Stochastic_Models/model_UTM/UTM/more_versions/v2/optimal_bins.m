%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This function computes the optimal amount of bins for a histogram using
% the Freedman-Diaconis rule.
%
% Input:
% x   : vector or matrix with observed data.
% dim : (optional) returns the optimal number of bins over the dimension 
%       "dim" in a matrix by n x m.
%
% Output:
% nbin : optimal number of bins.
%
% Christian González - PUC. April/2016.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [nbin] = optimal_bins(x,dim)

[n m] = size(x);
N = [n m];
if min(n,m)==1
    matrix = 0;
else
    matrix = 1;
end


if nargin == 1
    if matrix ==1
        bin_size = 2*(quantile(x,0.75) - quantile(x,0.25)) / (n^(2/3));
        r = n;
    else
        bin_size = 2*(quantile(x,0.75) - quantile(x,0.25)) / (length(x)^(2/3));
        r = length(x);
    end
    nbin = round((max(x)-min(x))./bin_size);
elseif nargin == 2    
	bin_size = 2*(quantile(x,0.75,dim) - quantile(x,0.25,dim)) / (N(dim)^(2/3));
    nbin = round((max(x,[],dim)-min(x,[],dim))./bin_size);
    r = N(dim);
end

% Remove infinites
nbin(nbin==inf)=r;




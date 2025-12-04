%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% COMPUTES PARAMETERS FOR THE GAMMA PDF FROM A HISTOGRAM USING METHOD OF MOMENTS
%
% Input:  x: x values of the histogram.
%         N: vector with histogram values.
%
% Output: k     : shape parameter of gamma distribution.
%         lambda: scale parameter of gamma distribution.
%
% Note: the difference between this function and gamfit is that gamfit uses
% as input the sampled data. This function uses as input the histogram of
% the sampled data.
%
% Created by Christian González. PUC. August, 2015.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


function [k lambda] = gamma_parameters_from_histogram(x,N)

mu = sum(x.*N)/sum(N);                    % mean
S2 = sum(N.*(x-mu).^2)/(sum(N)-1);        % variance

k      = mu^2/S2;                         % parameter k (k=a from gamfit)
lambda = mu/S2;                           % parameter lambda (lambda=1/b from gamfit)

return


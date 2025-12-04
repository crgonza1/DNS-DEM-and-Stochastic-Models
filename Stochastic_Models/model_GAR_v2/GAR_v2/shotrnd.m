% GENERATE RANDOM NUMBERS DISTRIBUTING SHOT-NOISE
%
% Let a GAR(1) model:
%
%              X(t) = r*X(t-1) + eta
%
% where X distributes Gamma(K,beta):
%
%     Gamma_pdf = [(1/beta)^K * X^(K-1) *exp(-X/beta)]/gamma(K)
%
% 
% INPUTS
% beta           : scale parameter of gamma distribution
% K              : shape parameter of gamma distribution
% r              : autoregression coefficient of GAR(1) model
% dim (optional) : size of the output
%
% OUTPUT
% eta            : random number
%
% E.g. eta = shotrnd(K,beta,r,[M,N]);
%      eta = shotrnd(10,5,0.5);
%      eta = shotrnd(10,5,0.5,[3,2]);
%
% Christian González M. - PUC. Nov, 2016
% crgonza1@uc.cl
%
% Ref. Fernández and Salas (1990) - Gamma autoregressive models for
% stream-flow simulation
%


function eta = shotrnd(varargin)

% Assign variables
K    = varargin{1};
beta = varargin{2};
r    = varargin{3};

if nargin == 4
    dim = varargin{4};
else
    dim = [1,1];
end


% Generate random number following a shot-noise distribution
eta = zeros(dim);
M   = poissrnd(-K*log(r),dim);

for i=1:dim(1)
    for j=1:dim(2)
                
        if M(i,j) == 0
            eta(i,j) = 0;
        else

            U = rand(M(i,j),1);
            Y = exprnd(beta,[M(i,j),1]);            
            for n=1:M(i,j)               
               eta(i,j) = eta(i,j) + Y(n)*r^U(n);                 
            end    
        end
        
        
    end       
end

return
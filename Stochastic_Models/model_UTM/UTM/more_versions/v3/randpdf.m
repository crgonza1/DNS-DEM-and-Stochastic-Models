%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This function generates a random number distributing according to an
% input pdf.
%
% Input:
% X   : vector of observed values.
% pdf : vector with the pdf associated to X.
% a   : number of rows of the output vector (optional).
% b   : number of columns of the output matrix (optional).
%
% Output:
% random_x : random value with the same distribution than the input pdf.
%
% Christian González - PUC. April/2016.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [random_x] = randpdf(X,pdf,a,b)

% Clean the zero-probability points
c=1;
for i=1:length(X)
    if pdf(i)>0
        pdf_aux(c) = pdf(i);
        X_aux(c)   = X(i);
        c = c+1;
    end
end
pdf = pdf_aux;
X   = X_aux;

% Create the pdf-space
cum(1) = pdf(1);
for i=2:length(pdf)
   cum(i) = pdf(i) + cum(i-1);
end
cum = [0 cum];

% Generate a random number
if nargin == 2
    r = cum(end)*rand;
    a = 1;
    b = 1;
elseif nargin == 3
    r = cum(end)*rand(a);
    b = 1;
elseif nargin == 4
    r = cum(end)*rand(a,b);
end

% Loop over the output random matrix (a x b)
for j=1:a
    for k=1:b
        % Find the index where cum(i) <= r <cum(i+1)
        i = find(cum<=r(j,k));  i = i(end);
        if i==size(cum)
            i = i-1;
        end
        
        % Assign the random value which distributes as the input pdf
        random_x(j,k) = X(i);
    end
end
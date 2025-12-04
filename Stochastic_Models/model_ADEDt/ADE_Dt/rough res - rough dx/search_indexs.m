%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%  Search the indexs a and b that satisfy V(a)<= x < V(b).
%
% Inputs: 
% V : sorted vector (from the smallest to the largest value).
% x : number to evaluate the indexs in the vector.
%
% Outputs:
% a : left index of the vector.
% b : right index of the vector.
%
% E.g. 
% >> V = [2 4 6 8];
% >> x = 4.1;
% >> [a,b] = search_indexs(V,x)
% a = 
%     2
% b =
%     3
%
% Christian González- PUC. May-2016.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [a,b] = search_indexs(V,x)

if x<V(1) | x>V(end)
   a = NaN;
   b = NaN;
   disp('Warning: the searched number is outside of the vector')
elseif x==V(1)
   a = 1;
   b = 2;   
elseif x==V(end)
   a = length(V)-1;
   b = length(V);
else
    k = find(V<=x); k=max(k);
    a = k;
    b = k+1;
end

return


% Remove values smaller than a given threshold.

function [x y] = remove_near_zero(x,y,threshold)

if length(x)~=1
    dx = mean(diff(x));
    
    % Remove small values
    x(abs(y)<threshold) = [];
    y(abs(y)<threshold) = [];
    
    % Generate new grid
    x_aux = min(x):dx:max(x);
    y = interp1(x,y,x_aux);
    x = x_aux;
    
    % Remove vestigials NaN
    x(isnan(y))=[];
    y(isnan(y))=[];
end

return
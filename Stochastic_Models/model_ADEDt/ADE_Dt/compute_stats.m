function [stats] = compute_stats(C,x,stats)

lt = size(C,1);

% Preallocate variables
xp_avg_ADE  = zeros(lt,1);
xp_var_ADE  = zeros(lt,1);
xp_skew_ADE = zeros(lt,1);
xp_kur_ADE  = zeros(lt,1);

% Compute statistics
for i=1:lt
    xp_avg_ADE(i,1)  = trapz(x,(C(i,:).*x))/trapz(x,C(i,:));
    xp_var_ADE(i,1)  = trapz(x,C(i,:).*(x - xp_avg_ADE(i)).^2)/trapz(x,C(i,:));
    xp_skew_ADE(i,1) = (trapz(x,(C(i,:).*x.^3))/trapz(x,C(i,:)) - 3*xp_avg_ADE(i)*xp_var_ADE(i) - xp_avg_ADE(i)^3)/sqrt(xp_var_ADE(i))^3;    
    xp_kur_ADE(i,1)  = (trapz(x,(C(i,:).*x.^4))/trapz(x,C(i,:)) - 4*xp_avg_ADE(i)*trapz(x,(C(i,:).*x.^3))/trapz(x,C(i,:)) + 6*xp_avg_ADE(i)^2*xp_var_ADE(i) + 3*xp_avg_ADE(i)^4)/xp_var_ADE(i)^2;
end

stats = [stats; [xp_avg_ADE xp_var_ADE xp_skew_ADE xp_kur_ADE]];

return

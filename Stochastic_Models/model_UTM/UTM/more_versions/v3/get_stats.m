function [t,xp_mean,xp_var,xp_skew,xp_kur] = get_stats(np,max_sim_time,dt,xi,pdf_xi)

% Define the time discretization
t = 0:dt:max_sim_time;

% Generate random length steps for each dt
rand_xi = randpdf(xi,pdf_xi,np,length(t)-1);

% Compute the particles position for each particle
aux       = cumsum(rand_xi,2);
positions = [zeros(np,1) aux];

% Compute the particles statistics as a function of the time
xp_mean = mean(positions);
xp_var  = var(positions,0,1);
xp_skew = skewness(positions,1,1);
xp_kur  = kurtosis(positions,1,1);

return
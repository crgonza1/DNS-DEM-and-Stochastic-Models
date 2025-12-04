function [U] = get_U(t)

global global_time_sim
global global_mean_vp_sim

U = interp1(global_time_sim,global_mean_vp_sim,t);

return

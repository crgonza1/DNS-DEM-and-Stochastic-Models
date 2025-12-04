function [D] = get_D(t)

global global_time_sim
global global_D_coef

D = interp1(global_time_sim,global_D_coef,t);

return

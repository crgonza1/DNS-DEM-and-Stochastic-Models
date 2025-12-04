function [c,f,s] = define_eq(x,t,C,DCDx)

D = get_D(t);
U = get_U(t);

%D = 0.0002;% +0.002*t;  % diffusion coefficient
%U = 0;%0.005;%0.01*t;       % mean velocity (dmean_xp_sim/dt)              


c = 1;
f = D*DCDx;
s = -U*DCDx;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Solve Diffusion - Advection Equation
% Input : x and t
% Output: C
%
% See more info:
% https://www.mathworks.com/help/matlab/math/partial-differential-equations.html?requestedDomain=www.mathworks.com#brfhed9-1
% Example: "edit pdex1.m"
%
% Christian González, April 2016.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [C] = ADE(varargin)

t = varargin{1};
x = varargin{2};

if nargin==2
    IC = @define_init_cond;
elseif nargin>=3
	IC = varargin{3};
end

% Concentration: C(t,x)
sol = pdepe(0,@define_eq,IC,@define_bound_cond,x,t);
C   = sol(:,:,1);

return
































return






% A surface plot is often a good way to study a solution.
figure;
surf(x,t,C);    
title('Numerical solution computed with 20 mesh points.');
xlabel('Distance x');
ylabel('Time t');


% A solution profile can also be illuminating.
figure;
plot(x,C(end,:),'o');
title('Solutions at t = 2.');
legend('Numerical, 20 mesh points','Analytical',0);
xlabel('Distance x');
ylabel('C(x,2)');
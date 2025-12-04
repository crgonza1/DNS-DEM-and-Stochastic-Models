function [t,xp_mean,xp_var,xp_skew,xp_kur,BTC] = model_AR(tmax,dt,K,rho,xi,pdf_xi,diameter,up,hist_up)

% Set parameters:
NP     = 1000000;                                                          % Number of particles
t      = [0:dt:tmax]';
s      = zeros(NP,1);
BTC    = zeros(length(t),5);
timer1 = cputime;
method = 1;                                                                % 1: use xi-pdf_xi data for computing mu_s.
                                                                           % 2: use up-hist data for computing mu_s.

% Compute the correlation time scale and the lag-dt autocorrelation
r = rho(2);


% Define the mean and standard deviation of the noise term of the particles steps
if method == 1
    mu_s    = trapz(xi,xi.*pdf_xi)/trapz(xi,pdf_xi);
    sigma_s = sqrt(trapz(xi,pdf_xi.*(xi-mu_s).^2)/trapz(xi,pdf_xi));
    
elseif method == 2
    mean_up  = trapz(up,up.*hist_up)/trapz(up,hist_up);
    var_up   = trapz(up,hist_up.*(up-mean_up).^2)/trapz(up,hist_up);
    mu_s     = mean_up*dt;
    sigma_s  = sqrt(var_up)*dt;
end


% Preallocate variables
xp_mean = zeros(length(t),1);
xp_var  = zeros(length(t),1);
xp_skew = zeros(length(t),1);
xp_kur  = zeros(length(t),1);


% Compute the statistics and BTCs
for i=1:length(t)
    
    % Definition of the noise
    eta = normrnd(mu_s*(1-r),sigma_s*sqrt(1-r^2),[NP,1]);
    
    % Sum of the steps
    if i==1
        xp = zeros(NP,1);
    elseif i==2
        xp    = normrnd(mu_s,sigma_s,[NP,1]);
        s_old = xp;
    else
        s  = r*s_old + eta;
        xp = xp + s;
        s_old = s;
    end
    
    % Compute the stats
    xp_mean(i,1) = mean(xp);
    xp_var(i,1)  = var(xp);
    xp_skew(i,1) = skewness(xp);
    xp_kur(i,1)  = kurtosis(xp);    
    
    % Add a new point to the BTCs
    for j=1:5
        dxups = (10^(j-2))*diameter;
        if min(xp)<=dxups && dxups <=max(xp)
            [N X]    = hist(xp,100);
            
            % Soft curve
            window = 1;
            N_aux = N;
            for k=1+window:length(N)-window
                N_aux(k) = mean(N(k-window:k+window));    
            end
            N = N_aux;
                        
            BTC(i,j) = interp1(X,N/trapz(X,N), dxups);
        end
    end
        
    % Display message
    if mod(i,floor(length(t)/100))==0
        timer2 = cputime;
        disp([num2str(round(100*i/length(t))),'% completed (current proccess took ',num2str(timer2-timer1),' s).'])
        timer1 = timer2;
    end
    
end

return
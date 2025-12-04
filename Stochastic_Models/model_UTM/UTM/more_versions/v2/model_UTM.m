function [t,xp_mean,xp_var,xp_skew,xp_kur,btc1,btc2,btc3,btc4,btc5] = model_UTM(np,max_sim_time,dt,xi,pdf_xi,diameter)

% Define the time discretization
t = 0:dt:max_sim_time;

% Generate random length steps for each dt
rand_xi      = randpdf(xi,pdf_xi,np,length(t)-1);

% Compute the particles position for each particle
aux          = cumsum(rand_xi,2);
positions    = [zeros(np,1) aux];

% Compute the particles statistics as a function of the time
xp_mean = mean(positions);
xp_var  = var(positions,0,1);
xp_skew = skewness(positions,1,1);
xp_kur  = kurtosis(positions,1,1);


% Get arrival travel times (for computing BTCs)
dxups = diameter*[0.1 1 10 100 1000];
btc_data = nan(np,length(dxups));
for i=1:np
           
    pos_org = positions(i,:);
    pos_aux = positions(i,:);
    t_aux   = t;
    
    % Clean the backward movements
    for j=2:length(t)
        if pos_org(j)<=pos_org(j-1)
            pos_org(j) = pos_org(j-1);
        end
    end
       
    % Find the time when the traveled distance is equal to dxups(k)
    for k=1:length(dxups)
        if max(pos_aux)>=dxups(k)
            x0            = intersections(t,pos_org,t,dxups(k)*ones(size(t)));
            btc_data(i,k) = min(x0);
        end  
    end
end

% Compute BTCs
for i=1:length(dxups)
    data_aux                  = btc_data(:,i);
    data_aux(isnan(data_aux)) = [];                                        % Remove NaN from the traveling times
    nbin                      = optimal_bins(data_aux);                    % Compute the optimal amount of bins of the BTCs histograms
    [N tbtc]                  = hist(data_aux,nbin);                       % Compute the histogram of the BTCs
    BTC{i}                    = [tbtc' N'/trapz(tbtc,N)];    
    if isempty(BTC{i})
        BTC{i} = [NaN NaN];
    end    
end

% Assigment of the BTCs
btc1 = BTC{1};
btc2 = BTC{2};
btc3 = BTC{3};
btc4 = BTC{4};
btc5 = BTC{5};

return
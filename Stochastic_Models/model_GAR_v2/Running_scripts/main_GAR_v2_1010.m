clear all
close all
clc


% Inputs ******************************************************************
diameter     = 0.0008;
path_inputs  = '/scratch365/cgonza13/Paper_II/rhop_1010';
path_outputs = '/scratch365/cgonza13/Upscaling_models/Results/GAR_v2/rhop_1010';
path_model   = '/scratch365/cgonza13/Upscaling_models/model_GAR_v2/GAR_v2';
direction    = 'x';
upmean       = 0.033976028095302;
num_part     = 100000;                                                      % Number of particles per processor
nproc        = 16;


% Load data ***************************************************************
disp('Loading data...')
data_hist_tvelx = load([path_inputs,'/pasted_histograms/tvelx_1_histogram.txt']);
data_dispersion = load([path_inputs,'/part_dispersion.txt']);
data_up         = load([path_inputs,'/pasted_histograms/up_histogram.txt']);
data_corr       = load([path_inputs,'/autocorr_tvelx.txt']);


% Assign variable *********************************************************
time_end = data_dispersion(end,1);
dt       = data_dispersion(2,1)-data_dispersion(1,1);
dt_ups   = max(1,floor((0.1*diameter/upmean) /dt)) * dt;
xi       = data_hist_tvelx(:,1)*dt_ups;
pdf_xi   = data_hist_tvelx(:,2)/trapz(xi,data_hist_tvelx(:,2));
up       = data_up(:,1);
hist_up  = data_up(:,2);
K        = data_corr(:,1);
rho      = data_corr(:,2);


% Check positivity of the total displacement
xi = check_xi(xi,pdf_xi,(data_dispersion(end,2)/time_end)*dt_ups);

clear data_hist_tvelx data_dispersion


% Call function ***********************************************************
disp('Calling GAR_v2 model...')
addpath(path_model)
tic


% Single processor computing
if nproc == 1
    [t,xp_mean,xp_var,xp_skew,xp_kur,BTC] = model_GAR_v2(time_end,dt_ups,K,rho,xi,pdf_xi,diameter,up,hist_up,num_part);

else
   % Parallel computing
    delete(gcp('nocreate'))
    parpool('local', nproc);
    
    % Parallel loop
    parfor myid = 1:nproc
        [t(:,myid),xp_mean(:,myid),xp_var(:,myid),xp_skew(:,myid),xp_kur(:,myid),BTC(:,:,myid)] = model_GAR_v2(time_end,dt_ups,K,rho,xi,pdf_xi,diameter,up,hist_up,num_part);
    end
    
    % Merge data from the different processors
    t       = t(:,1);
    xp_mean = mean(xp_mean,2);
    xp_var  = mean(xp_var,2);
    xp_skew = mean(xp_skew,2);
    xp_kur  = mean(xp_kur,2);
    BTC     = mean(BTC,3);
    
    delete(gcp('nocreate'))
end

disp(['The model finished in ',num2str(toc),' s.'])
rmpath(path_model)


% Remove NaN
xp_mean(isnan(xp_mean)) = 0;
xp_var(isnan(xp_var))   = 0;
xp_skew(isnan(xp_skew)) = 0;
xp_kur(isnan(xp_kur))   = 0;
BTC(isnan(BTC))         = 0;


% Write outputs ***********************************************************
disp('Writing outputs...')
fileID = fopen([path_outputs,'/',direction,'/','stats.plt'],'w');
fprintf(fileID,'%12.8e %12.8e %12.8e %12.8e %12.8e\n',[t(:) xp_mean(:) xp_var(:) xp_skew(:) xp_kur(:)]');
fclose(fileID);

fileID = fopen([path_outputs,'/',direction,'/','btc1.plt'],'w');
fprintf(fileID,'%12.8e %12.8e\n',[t(:) BTC(:,1)]');
fclose(fileID);

fileID = fopen([path_outputs,'/',direction,'/','btc2.plt'],'w');
fprintf(fileID,'%12.8e %12.8e\n',[t(:) BTC(:,2)]');
fclose(fileID);

fileID = fopen([path_outputs,'/',direction,'/','btc3.plt'],'w');
fprintf(fileID,'%12.8e %12.8e\n',[t(:) BTC(:,3)]');
fclose(fileID);

fileID = fopen([path_outputs,'/',direction,'/','btc4.plt'],'w');
fprintf(fileID,'%12.8e %12.8e\n',[t(:) BTC(:,4)]');
fclose(fileID);

fileID = fopen([path_outputs,'/',direction,'/','btc5.plt'],'w');
fprintf(fileID,'%12.8e %12.8e\n',[t(:) BTC(:,5)]');
fclose(fileID);


return
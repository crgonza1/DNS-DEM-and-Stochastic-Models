clear all
close all
clc


% Inputs ******************************************************************
diameter     = 0.0008;
path_inputs  = '/scratch365/cgonza13/Paper_II/rhop_1347';
path_outputs = '/scratch365/cgonza13/Upscaling_models/Results/UTM/rhop_1347';
path_model   = '/scratch365/cgonza13/Upscaling_models/model_UTM/UTM';
direction    = 'x';
upmean       = 0.018442584075884;


% Load data ***************************************************************
disp('Loading data...')
data_hist_tvelx = load([path_inputs,'/pasted_histograms/tvelx_1_histogram.txt']);
data_dispersion = load([path_inputs,'/part_dispersion.txt']);


% Assign variable *********************************************************
time_end = data_dispersion(end,1);
dt       = data_dispersion(2,1)-data_dispersion(1,1);
dt_ups   = max(1,floor((0.1*diameter/upmean) /dt)) * dt;
xi       = data_hist_tvelx(:,1)*dt_ups;
pdf_xi   = data_hist_tvelx(:,2)/trapz(xi,data_hist_tvelx(:,2));


% Check positivity of the total displacement
xi = check_xi(xi,pdf_xi,(data_dispersion(end,2)/time_end)*dt_ups);

clear data_hist_tvelx data_dispersion


% Call function ***********************************************************
disp('Calling UTM model...')
addpath(path_model)
tic
[t,xp_mean,xp_var,xp_skew,xp_kur,BTC] = model_UTM(time_end,dt_ups,xi,pdf_xi,diameter);
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
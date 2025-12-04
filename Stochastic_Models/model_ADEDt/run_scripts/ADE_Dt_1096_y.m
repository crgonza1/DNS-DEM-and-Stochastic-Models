clear all
close all
clc


% Inputs ******************************************************************
diameter     = 0.0008;
direction    = 'y';
path_inputs  = '/scratch365/cgonza13/Paper_II/rhop_1096';
path_outputs = '/scratch365/cgonza13/Upscaling_models/Results/ADE_Dt/rhop_1096';
path_model   = '/scratch365/cgonza13/Upscaling_models/model_ADEDt/ADE_Dt';
dx           = 0.1*diameter;
% *************************************************************************


% Load data ***************************************************************
disp('Loading data...')
data_dispersion = load ([path_inputs,'/part_dispersion.txt']);


% Assign variable *********************************************************
time        = data_dispersion(:,1);
if strcmp(direction,'x')
    meanpos = data_dispersion(:,2);
    varpos  = data_dispersion(:,4);
elseif strcmp(direction,'y')
    meanpos = data_dispersion(:,3);
    varpos  = data_dispersion(:,5);    
else
    disp('Please specify a valid direction.')
    return
end


% Call function ***********************************************************
disp('Calling ADE D(t) model...')
addpath(path_model)
tic
[time_ADE,meanpos_ADE,varpos_ADE,skew_ADE,kur_ADE,BTCs,D_coef] = model_ADE_Dt(time,meanpos,varpos,diameter,dx);
disp(['The model finished in ',num2str(toc),' s.'])
rmpath(path_model)


% Write outputs ***********************************************************
disp('Writing outputs...')
fileID = fopen([path_outputs,'/',direction,'/','stats.txt'],'w');
fprintf(fileID,'%12.8e %12.8e %12.8e %12.8e %12.8e\n',[time_ADE meanpos_ADE varpos_ADE skew_ADE kur_ADE]');
fclose(fileID);

fileID = fopen([path_outputs,'/',direction,'/','Dcoef.txt'],'w');
fprintf(fileID,'%12.8e %12.8e\n',[time_ADE D_coef]');
fclose(fileID);

fileID = fopen([path_outputs,'/',direction,'/','BTCs.txt'],'w');
fprintf(fileID,'%12.8e %12.8e %12.8e %12.8e %12.8e %12.8e\n',[time_ADE BTCs(:,1) BTCs(:,2) BTCs(:,3) BTCs(:,4) BTCs(:,5)]');
fclose(fileID);

%save([path_outputs,'/',direction,'/','C.mat'],'C')
%save([path_outputs,'/',direction,'/','x_axis.mat'],'x_ADE')
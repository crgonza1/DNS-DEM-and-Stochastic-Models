clear all
close all
clc


% Inputs ******************************************************************
diameter     = 0.0008;
direction    = 'x';
path_inputs  = '/scratch365/cgonza13/Paper_II/rhop_1347';
path_outputs = '/scratch365/cgonza13/Upscaling_models/Results/ADE_Dt/rhop_1347';
path_model   = '/scratch365/cgonza13/Upscaling_models/model_ADEDt/ADE_Dt';
dx           = 0.1*diameter;
% *************************************************************************



% BUG *********************************************************************
path_inputs  = '/Volumes/Elements/Khript/Itzamna/Paper_II/rhop_1080';
path_outputs = '/Users/Christian/Desktop/Khript/Google_Drive/Investigaciones/upscaling_models/src';
path_model   = '/Users/Christian/Desktop/Khript/Google_Drive/Investigaciones/upscaling_models/DEFINITIVOS/model_ADEDt/ADE_Dt';
% BUG *********************************************************************



% Load data ***************************************************************
disp('Loading data...')
data_dispersion = load ([path_inputs,'/part_dispersion.txt']);

data_dispersion = data_dispersion(1:50,:);% BUG


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


% BUG ******************
figure
subplot(2,2,1)
plot(data_dispersion(:,1),data_dispersion(:,2),'.',time_ADE,meanpos_ADE,'r')
subplot(2,2,2)
plot(data_dispersion(:,1),data_dispersion(:,4),'.',time_ADE,varpos_ADE,'r')
subplot(2,2,3)
plot(data_dispersion(:,1),data_dispersion(:,6),'.',time_ADE,skew_ADE,'r')
subplot(2,2,4)
plot(data_dispersion(:,1),data_dispersion(:,8),'.',time_ADE,kur_ADE,'r')

figure
subplot(2,3,1)
plot(time_ADE,BTCs(:,1))
subplot(2,3,2)
plot(time_ADE,BTCs(:,2))
subplot(2,3,3)
plot(time_ADE,BTCs(:,3))
subplot(2,3,4)
plot(time_ADE,BTCs(:,4))
subplot(2,3,5)
plot(time_ADE,BTCs(:,5))
% BUG ******************


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

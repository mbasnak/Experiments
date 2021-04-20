%code to check phase and stim convention

clear all; close all;

parentDir = 'Z:\Wilson Lab\Mel\Experiments\Uncertainty\Exp25\data\Experimental\two_ND_filters_3_contrasts\20201015_60D05_7f';
sid = 24;
tid = 0;

% find the ROI_analysis folder
roi_analysis_dir = [parentDir, '\2p\ROI_analysis\'];

% find the behavior ball data folder
ball_dir = [parentDir, '\ball\'];

% create the data analysis folder
data_analysis_dir = [parentDir, '\analysis\'];
if(~exist(data_analysis_dir, 'dir'))
    mkdir(data_analysis_dir);
end

%% Get ROI_analysis file, load

expression = ['*sid_' num2str(sid) '_tid_' num2str(tid) '*.mat'];
roi_analysis_file = dir(fullfile(roi_analysis_dir, expression));
load(fullfile(roi_analysis_dir, roi_analysis_file.name)); %load ROI data (dff, zscore,...)

%Sort by proper ROI order using PB glomeruli order (this way the 1st
%will be the glomerulus on the far right, which will get plotted on the
%top of the imagesc)
order = [8,9,7,10,6,11,5,12,4,13,3,14,2,15,1,16];

for i = 1:length(order)
    roi_data(i).order = order(i);%add field with the ROI order
end
T = struct2table(roi_data); % convert the struct array to a table
sortedT = sortrows(T, 'order'); % sort the table by ROI order
roi_data = table2struct(sortedT);

%% Import behavior file

ball_file = dir(fullfile(ball_dir, expression));
ballData = load(fullfile(ball_dir, ball_file.name)); %load ballData
runobjFile = dir(fullfile([ball_dir,'\runobj\'], ['*_sid_' num2str(sid) '_*']));
load(fullfile([ball_dir '\runobj\'], runobjFile.name)); %load run_obj


%% Create dF/F matrix (1-16)

volumes = length(roi_data(1).dff); %get the volume number

dff_matrix = zeros(16, volumes);
for i = 1:16
    dff_matrix(i, :) = roi_data(i).dff; %I changed this because the field 'num' does not exist - MB 20190917
end

figure,
subplot(4,1,1)
imagesc(dff_matrix)
colormap(gray)
yticks(1:2:16);
yticklabels({'8R','6R','4R','2R','1L','3L','5L','7L'});
ylabel('PB glomerulus');
title('EPG activity in the PB');

%% Calculate power spectrum and phase using Fourier transform

transform = fft(dff_matrix, 16, 1); 
transform(1,:) = [];

n = 16;
phase = angle(transform(1:8,:));  %take the phase of the FFT of period 8 glomeruli
phase = squeeze(phase(2,:)); 

%add phase to plot
subplot(4,1,2)
plot(wrapTo360(rad2deg(phase)))
ylim([0 360]);
xlim([0,length(phase)]);
title('Phase obtained from Fourier transform');

s = .5.^(1:.025:5); % create a geometric series for the periodogram -- this is the FREQUENCY
pxx = periodogram(dff_matrix, [], s, 1, 'power'); % POWER SPECTRA
position_eight = find(s == .125); % Look for Period of 8 glomeruli
power_value = pxx(position_eight, :);

%% Calculate PVA

%Separate the PB into the complimentary halves
left_dff = dff_matrix([9,16,15,14,13,12,11,10],:);
right_dff = dff_matrix([8,7,6,5,4,3,2,1],:);

%Average both halves
mean_dff = (left_dff + right_dff)./2;

%Shift the mean to obtain the top of the EB (i.e., glomeruli 5 of the PB))
%at the top
mean_dff_EB = circshift(mean_dff,4);
dff_pva = circ_mean(repmat([pi/8:pi/4:15*pi/8], size(mean_dff_EB,2),1), mean_dff_EB', 2);

subplot(4,1,3)
imagesc(flip(mean_dff_EB))
colormap(gray)
yticks(1:2:8);
yticklabels({'8','6','4','2'});
ylabel('EB tiles');
title('EPG activity transformed to EB coordinates');

subplot(4,1,4)
plot(wrapTo360(rad2deg(dff_pva)));
xlim([0 length(dff_pva)]);
ylim([0 360]);
title('Phase as PVA');


%% Add behavior

figure,
plot(360-ballData.trial_bdata(:,5)*36)

clear all; close all;
%% Load data

load('Z:\Wilson Lab\Mel\Experiments\Uncertainty\Exp30\data\20210312_60D05_7f_fly3\analysis\analysis_sid_0_tid_0.mat')
%load analysis_sid_0_tid_0.mat

%% Method 1: old method, without shifting the data, and computing bump magnitude as max-min

figure('Position',[100 100 1400 600]),
subplot(3,1,1)
imagesc(data.dff_matrix)
colormap(flipud(gray))
title('EPG activity');

%compute phase
transform = fft(data.dff_matrix, 16, 1);
transform(1,:) = []; 
phase = angle(transform(1:8,:));  
phase_value = squeeze(phase(2,:));
subplot(3,1,2)
plot(phase_value)
xlim([0 length(phase_value)]);
ylim([-pi pi]);
title('EPG phase');

%compute bump magnitude
bump_mag = max(data.mean_dff_EB)-min(data.mean_dff_EB);
subplot(3,1,3)
plot(bump_mag,'r')
xlim([0 length(phase_value)]);
title('BM max-min');
xlim([0 length(data.dff_matrix)]);


%% Method 2: using Anna's function with the shifted data according to corresponding PB halves

shifted_dff_matrix = data.dff_matrix([9,16,15,14,13,12,11,10,8,7,6,5,4,3,2,1],:);

figure('Position',[100 100 1400 600]),
subplot(3,1,1)
imagesc(data.dff_matrix)
colormap(flipud(gray))
title('EPG activity');

%compute phase and magnitude
[phase,amplitude] = FT_Anna(shifted_dff_matrix);
phase_value_A = phase(3,:);
bump_mag_A = amplitude(3,:);

subplot(3,1,2)
plot(phase_value_A)
xlim([0 length(phase_value_A)]);
ylim([-pi pi]);
title('EPG phase');

%compute bump magnitude
subplot(3,1,3)
plot(bump_mag_A,'r')
xlim([0 length(phase_value_A)]);
title('BM max-min');
xlim([0 length(data.dff_matrix)]);



%% Method 3: using new method with the shifted data according to corresponding PB halves

figure('Position',[100 100 1400 600]),
subplot(3,1,1)
imagesc(data.dff_matrix)
colormap(flipud(gray))
title('EPG activity');

%compute phase
transform = fft(shifted_dff_matrix, 16, 1);
transform(1,:) = []; 
phase = angle(transform(1:8,:));  
phase_value_3 = squeeze(-phase(2,:));
subplot(3,1,2)
plot(phase_value_3)
xlim([0 length(phase_value_3)]);
ylim([-pi pi]);
title('EPG phase');

%compute bump magnitude
amplitude =  2*abs(transform/16);
bump_mag_3 = amplitude(2,:);
subplot(3,1,3)
plot(bump_mag_3,'r')
xlim([0 length(phase_value_3)]);
title('BM max-min');
xlim([0 length(data.dff_matrix)]);


%% Phase value comparison

figure,
subplot(3,1,1)
plot(phase_value)
title('Phase value with previous method');
xlim([0 length(data.bump_magnitude)]);

subplot(3,1,2)
plot(phase_value_A)
title('Phase value with Anna method');
xlim([0 length(data.bump_magnitude)]);

subplot(3,1,3)
plot(phase_value_3)
title('Phase value with current method');
xlim([0 length(data.bump_magnitude)]);


%% Bump magnitude comparison

figure,
subplot(3,1,1)
plot(bump_mag)
title('BM with previous method');
xlim([0 length(data.bump_magnitude)]);

subplot(3,1,2)
plot(bump_mag_A)
title('BM with Anna method');
xlim([0 length(data.bump_magnitude)]);

subplot(3,1,3)
plot(bump_mag_3)
title('BM with current method');
xlim([0 length(data.bump_magnitude)]);

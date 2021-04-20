%compute magnitude of first fft component
clear 
load analysis_sid_1_tid_0.mat
signal = data.dff_matrix;
tic
transform = fft(signal, 16, 1); %apply a fast fourier transform to the dff data   
transform(1,:) = []; % remove 0th order (actually not needed, but ok if done)
phase_fft = angle(transform(1:8,:)); % compute phase of components 1 - 8
phase_value = squeeze(-phase_fft(2,:)); % the phase of the !2nd! component - what about shift?
toc
tic
[phase,amplitude] = FT_Anna(signal);
toc

%plot things
figure(1), clf
%heatmap of EPG activity
subplot(3,1,1)
imagesc(data.dff_matrix)
colormap(flipud(gray))
title('EPG activity');

%EPG activity phase
subplot(3,1,2)
plot(phase_value)
hold on
plot(-phase(3,:))
xlim([0 length(transform)])
title('Phase value');

%magnitude of the first component
subplot(3,1,3)
plot(abs(transform(2,:)./(16/2) )) % needed correction
hold on
plot(amplitude(3,:))
xlim([0 length(transform)])
title('Value of second Fourier component');


% ANNA ADDED: previously, you used real(), which does not give you the amplitude, but
% only the real part of the imaginary number.
% % %absolute magnitude of the first component
% % subplot(4,1,4)
% % plot(abs(real(transform(1,:)))) 
% % xlim([0 length(transform)])
% % title('Magnitude of first Fourier component');

%I think this data does not look like the bump magnitude at all, so I'm
%next repeating this analysis without removing that first original
%component

%% Averaged data, no shift

signal = data.mean_dff_EB;
[phase_mean,amplitude_mean] = FT_Anna(signal);

%plot things
figure(2), clf
%heatmap of averaged EPG activity
subplot(3,1,1)
imagesc(signal)
colormap(flipud(gray))
title('EPG activity, mean across glomeruli');

%EPG activity phase
subplot(3,1,2)
plot(phase_value) % as raw data as comparin
hold on
plot(phase_mean(2,:))
xlim([0 length(transform)])
title('Phase value')
legend('raw','mean')

%magnitude of the first component
subplot(3,1,3)
plot(abs(transform(2,:)./(16/2) )) % needed correction
hold on
plot(amplitude_mean(2,:))
xlim([0 length(transform)])
title('Value of first Fourier component')
legend('raw','mean')

%% %% Only data from glomeroli 1-8

signal = data.dff_matrix(1:8,:);
[phase_half,amplitude_half] = FT_Anna(signal);

%plot things
figure(3), clf
%heatmap of averaged EPG activity
subplot(3,1,1)
imagesc(signal)
colormap(flipud(gray))
title('EPG activity, only glomeroli 1 - 8');

%EPG activity phase
subplot(3,1,2)
plot(phase_value) % as raw data as comparin
hold on
plot(-phase_half(2,:))
xlim([0 length(transform)])
title('Phase value')
legend('raw','half')

%magnitude of the first component
subplot(3,1,3)
plot(abs(transform(2,:)./(16/2) )) % needed correction
hold on
plot(amplitude_half(2,:))
xlim([0 length(transform)])
title('Value of first Fourier component')
legend('raw','half')

%% 
max(mod(abs(amplitude(3,:)-abs(transform(2,:)./(16/2))), 2 * pi))
max(mod(abs(amplitude_mean(2,:)-abs(transform(2,:)./(16/2))), 2 * pi))
max(mod(abs(amplitude_half(2,:)-abs(transform(2,:)./(16/2))), 2 * pi))
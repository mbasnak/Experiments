%compute magnitude of first fft component
transform = fft(data.dff_matrix, 16, 1); %apply a fast fourier transform to the dff data   
transform(1,:) = []; %Jenny usually removes the first component before computing the phase value
phase = angle(transform(1:8,:));
phase_value = squeeze(phase(2,:));

%plot things
figure,
%heatmap of EPG activity
subplot(4,1,1)
imagesc(data.dff_matrix)
colormap(flipud(gray))
title('EPG activity');

%EPG activity phase
subplot(4,1,2)
plot(phase_value)
xlim([0 length(transform)])
title('Phase value');

%magnitude of the first component
subplot(4,1,3)
plot(real(transform(1,:)))
xlim([0 length(transform)])
title('Value of first Fourier component');

%absolute magnitude of the first component
subplot(4,1,4)
plot(abs(real(transform(1,:))))
xlim([0 length(transform)])
title('Magnitude of first Fourier component');

%I think this data does not look like the bump magnitude at all, so I'm
%next repeating this analysis without removing that first original
%component

%% if I don't remove that first component

transform2 = fft(data.dff_matrix, 16, 1);  
phase2 = angle(transform2(1:8,:));
phase_value2 = squeeze(phase2(2,:));

%plot things
figure,
%heatmap of EPG activity
subplot(4,1,1)
imagesc(data.dff_matrix)
colormap(flipud(gray))
title('EPG activity');

%EPG activity phase
subplot(4,1,2)
plot(phase_value2)
xlim([0 length(transform)])
title('Phase value')

%magnitude of the first component
subplot(4,1,3)
plot(real(transform2(1,:)))
xlim([0 length(transform)])
title('Value of first Fourier component');

%absolute magnitude of the first component
subplot(4,1,4)
plot(abs(real(transform2(1,:))))
xlim([0 length(transform)])
title('Magnitude of first Fourier component');

%this looks more like it in terms of the fourier component magnitude

%% Overlaying the phase values

figure,
subplot(2,1,1)
imagesc(data.dff_matrix)
colormap(flipud(gray))
title('EPG activity');

subplot(2,1,2)
plot(phase_value)
hold on
plot(phase_value2)
xlim([0 length(transform)])
legend('Removing first component','Without removing it');
title('EPG phase values');

%the phase values don't look the same, and the one computed using Jenny's
%method is the one that looks correct!

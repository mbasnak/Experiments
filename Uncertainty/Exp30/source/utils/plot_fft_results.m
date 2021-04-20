
function plot_fft_results(dff_matrix,phase_value,amplitude,power_value)

figure
%epg activity
subplot(4,1,1)
imagesc(dff_matrix)
colormap(flipud(gray))
title('EPG activity');

%phase
subplot(4,1,2)
plot(phase_value,'k')
ylim([-pi pi]);
xlim([0 length(phase_value)]);
title('Phase value obtained from FFT');

%bump magnitude
subplot(4,1,3)
plot(amplitude,'r')
xlim([0 length(amplitude)]);
title('Bump magnitude obtained from FFT');

%ft power
subplot(4,1,4)
plot(power_value,'b')
xlim([0 length(power_value)]);
title('Power value of FFT');
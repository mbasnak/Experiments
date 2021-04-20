%code to get the data from the hdf5 file

clear all; close all;

path_name = 'Z:\Wilson Lab\Mel\Experiments\PathIntegration\Exp26\data\Debug\20200925_debug_fictrac\ball\'
file_name = [path_name,'hdf5_inbound_outbound_Closed_Loop_X_Closed_Loop_Y_20200925_130548_sid_8_tid_1000001'];


timestamp = h5read([file_name,'.hdf5'],'/time');
gain = h5read([file_name,'.hdf5'],'/gain');
voltage_out = h5read([file_name,'.hdf5'],'/output_voltage_x_gain');
panel_x = h5read([file_name,'.hdf5'],'/panel x');
intx = h5read([file_name,'.hdf5'],'/intx');
velx = h5read([file_name,'.hdf5'],'/velx');


%% Plot results

figure,
subplot(3,1,1)
plot(timestamp,velx)
xlim([0, timestamp(end)]);
ylim([min(velx),max(velx)]);
title('Change in the x vel from fictrac noise');

subplot(3,1,2)
plot(timestamp(velx>0),velx(velx>0))
xlim([0, timestamp(end)]);
ylim([min(velx),max(velx)]);
title('Keeping only positive changes in x vel');

subplot(3,1,3)
plot(timestamp,panel_x)
xlim([0, timestamp(end)]);
ylim([0 360]);
title('Change in the panels coming from fictrac drift');


%Make velx <0 be = 0 (which is what my code is doing by preventing the
%animals to go backwards

constrained_velx = velx;
constrained_velx(velx<0) = 0;

figure,
subplot(2,1,1)
plot(timestamp,constrained_velx)
title('Velx constrained to positive values');
xlim([0, timestamp(end)]);

subplot(2,1,2)
accumx = cumsum(constrained_velx);
plot(timestamp,accumx)
title('Drift in x position');
xlim([0, timestamp(end)]);

%Find the total and instantaneous drift
drift = accumx(end);
inst_drift = drift/length(timestamp);

%Plot the accumx position with a drift correction
for time = 1:length(timestamp)
    corrected_velx(time) = constrained_velx(time) - inst_drift;
end

figure, 
subplot(3,1,1)
plot(corrected_velx)

subplot(3,1,2)
plot(cumsum(corrected_velx))

%what the panel output woul look like in this case
inst_panel = corrected_velx*(360/2*pi);
subplot(3,1,3)
plot(wrapTo360(cumsum(inst_panel)));
ylim([0 360]);
title('Change in panels when correcting for the drift');
function [motor_pos] = process_motor_360(motor)

%code to process the data form the motor moving the closed-loop wind device

sampRate_new = 25; % sampling rate after downsampling (Hz)
settings = nidaq_settings;

downsampled.motor = resample_new(motor, sampRate_new, settings.sampRate);
motor_pos = downsampled.motor .* 2 .* pi ./ 10; % from voltage to radian

end
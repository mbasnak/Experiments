function settings = nidaq_settings

% Acquisition params
settings.sampRate = 4000;

% Processing settings
settings.cutoffFreq = 100;
settings.aiType = 'SingleEnded';
settings.sensorPollFreq = 50; 

settings.fictrac_x_DAQ_AI = 1;
settings.fictrac_y_DAQ_AI = 2; 
settings.fictrac_yaw_DAQ_AI = 3; % In this experiment, this is actually not yaw, but the panels y dimension
settings.fictrac_yaw_gain_DAQ_AI = 4;

settings.panels_x_DAQ_AI = 5;
settings.panels_y_DAQ_AI = 7;
settings.piezo = 8;
settings.motor_DAQ_AI = 9;
settings.wind_valve = 6;

settings.mfc_control = 10;
settings.mfc_monitor = 11;


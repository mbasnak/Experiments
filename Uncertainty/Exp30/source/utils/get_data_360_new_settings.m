function [smoothed, t, angle, flyPosRad, jump_frames, contrast] = get_data_360_new_settings(trial_time, trial_data)

%% Processes data from ball and panels for a 360 degree panel arena.

% DAQ channels and sensor settings, this will change if using a different
% setup
settings = sensor_settings;
settings.fictrac_yaw_gain_DAQ_AI = 1;
settings.fictrac_x_DAQ_AI = 5;
settings.fictrac_yaw_DAQ_AI = 3;
settings.fictrac_y_DAQ_AI = 2;
settings.panels_DAQ_AI = 4;
settings.panels_DAQ_AI_Y = 7;
settings.motor = 8;

% Assign fictrac data and compute position and velocities in proper units
data.ficTracIntx = trial_data(:,settings.fictrac_x_DAQ_AI);
data.ficTracAngularPosition_uncorrected = trial_data(:,settings.fictrac_yaw_DAQ_AI);
data.ficTracInty = trial_data(:,settings.fictrac_y_DAQ_AI);

% Assign fictrac data and compute position and velocities in proper units
data.ficTracIntx = trial_data(:,settings.fictrac_x_DAQ_AI);
data.ficTracAngularPosition_uncorrected = trial_data(:,settings.fictrac_yaw_DAQ_AI);
data.ficTracInty = trial_data(:,settings.fictrac_y_DAQ_AI);

%%  Correct the heading around the jumps (only for bar jump trials)

%1) Identify the changes in stim using the signal from the yPanels channel
panels_y = trial_data(:,settings.panels_DAQ_AI_Y);

%identify if the trial is a bar jump trial
if std(panels_y) > 2 
diff_panels_y = diff(panels_y);

%2) Find the jumps
h_to_l_120 = find(diff_panels_y>4 & diff_panels_y<6);
h_to_l_neg_120 = find(diff_panels_y>7);
l_to_l_120 = find(diff_panels_y<-4);
l_to_l_neg_120 = find(diff_panels_y>-4 & diff_panels_y<-2 & panels_y(1:end-1)>7);
l_to_h_0 = find(diff_panels_y>-4 & diff_panels_y<-2 & panels_y(1:end-1)<7);

%3) Plot to check
% figure,
% plot(panels_y)
% hold on
% plot(h_to_l_120,panels_y(h_to_l_120),'ro')
% plot(h_to_l_neg_120,panels_y(h_to_l_neg_120),'ko')
% plot(l_to_l_neg_120,panels_y(l_to_l_neg_120),'go')
% plot(l_to_l_120,panels_y(l_to_l_120),'co')
% plot(l_to_h_0,panels_y(l_to_h_0),'yo')

%4) Extract those that correspond to actual jumps in the stimulus position
positive_jump_frames = sort([h_to_l_120;l_to_l_120]);
negative_jump_frames = sort([h_to_l_neg_120;l_to_l_neg_120]);
jump_frames = sort([positive_jump_frames;negative_jump_frames]);
jump_voltages = [];
for jump = 1:length(jump_frames)
    if any(positive_jump_frames == jump_frames(jump))
        jump_voltages(jump) = -120*10/360; %this needs to be negative to match the plotting convention
    else
        jump_voltages(jump) = 120*10/360;
    end
end
cum_jump_voltages = cumsum(jump_voltages);

%The x position of the panels and the FicTrac output are "ignorant" of the
%bar jumps. The x position of the bar will move with the angular position
%of the fly, but the coordinate system changes every time the bar jumps.
%We need to make sure the xpos and heading of the fly are corrected to take
%this coordinate change into account.
panels_x_uncorrected = trial_data(:,settings.panels_DAQ_AI);
panels_x_corrected_uw = panels_x_uncorrected;
data.ficTracAngularPosition_uw = data.ficTracAngularPosition_uncorrected;

for i = 1:size(jump_frames)-1
    panels_x_corrected_uw(jump_frames(i)+1:jump_frames(i+1)) = panels_x_uncorrected(jump_frames(i)+1:jump_frames(i+1))+cum_jump_voltages(i);      
    data.ficTracAngularPosition_uw(jump_frames(i)+1:jump_frames(i+1)) = (data.ficTracAngularPosition_uncorrected(jump_frames(i)+1:jump_frames(i+1)))+cum_jump_voltages(i);
end
%Correct the data after the last jump
panels_x_corrected_uw(jump_frames(end)+1:end) = panels_x_uncorrected(jump_frames(end)+1:end)+cum_jump_voltages(end);
data.ficTracAngularPosition_uw(jump_frames(end)+1:end) = data.ficTracAngularPosition_uncorrected(jump_frames(end)+1:end)+cum_jump_voltages(end);

%wrap the data
panels_x_corrected = wrapTo10(panels_x_corrected_uw);
data.ficTracAngularPosition = wrapTo10(data.ficTracAngularPosition_uw);

% %Plot the original and corrected panel signals
% figure,
% plot(panels_x_uncorrected)
% for jump=1:length(jump_frames)
% xline(jump_frames(jump),'r','linewidth',3)
% end
% hold on
% plot(panels_x_corrected,'color',[0 0.6 0])
% 
% 
% %Check the jump size
% data.xPanelDeg = round((panels_x_corrected*360)/10); % Convert from what it reads in volts from the Ni-Daq to an X position in pixels in the panels
% for jump = 1:length(jump_frames)
%     figure,
%     plot(data.xPanelDeg(jump_frames(jump)-40000:jump_frames(jump)+40000))
%     hold on
%     jump_mag(jump) = rad2deg(circ_dist(deg2rad(data.xPanelDeg(jump_frames(jump)+1)),deg2rad(data.xPanelDeg(jump_frames(jump)-1))));
%     xline(40001,'r')
%     title(['Jump #',num2str(jump)]);
%     text(41000,200,num2str(jump_mag(jump)))
% end


% Store info on whether the contrast was high or low
contrast = [];
contrast(panels_y < 3) = 1; %high contrast
contrast(panels_y > 3) = 0;

else %if this is a normal trial
    panels_x_corrected = trial_data(:,settings.panels_DAQ_AI);
    data.ficTracAngularPosition = data.ficTracAngularPosition_uncorrected;
    jump_frames = [];
    contrast = [];
end
%% Compute velocity parameters
%%%%%%%%%%%%%%%

[smoothed] = singleTrialVelocityAnalysis9mm(data,settings.sampRate);
flyPosRad = smoothed.angularPosition;


%% 

% Get stimulus position
[angle] = process_panel_360(panels_x_corrected);

% Get time in proper units
[t] = process_time(trial_time);

end


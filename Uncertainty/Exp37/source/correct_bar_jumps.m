
function [smoothed, t, angle, flyPosRad, motor_pos] = correct_bar_jumps( trial_time, trial_data, num_frames_x)
%This functions takes the behavior data and corrects it taking into account
%the stimulus jumps (which it identifies using the panels y signal)

%% Set DAQ channels and sensor settings 

settings = sensor_settings;
settings.fictrac_yaw_gain_DAQ_AI = 4;
settings.fictrac_x_DAQ_AI = 1;
settings.fictrac_yaw_DAQ_AI = 3;
settings.fictrac_y_DAQ_AI = 2;
settings.panels_DAQ_AI = 5;
settings.panels_DAQ_AI_Y = 6;
settings.motor = 8;

ft_for = trial_data( :, settings.fictrac_x_DAQ_AI );
ft_yaw = trial_data( :, settings.fictrac_yaw_DAQ_AI );
ft_side = trial_data( :, settings.fictrac_y_DAQ_AI );
panels = trial_data( :, settings.panels_DAQ_AI );
panels_y_raw = trial_data( :, settings.panels_DAQ_AI_Y );
motor = trial_data( :, settings.motor );

%% Find and plot the jump frames

%Use the y dimension of the panels to find the jumps, taking its derivative
jumps = diff(panels_y_raw);
jumps(abs(jumps)>0.4 & abs(jumps)<1)=1;
jumps = round(jumps);
j = find(jumps); %indices of the actual bar jumps, taken from the y signal

% figure,
% suptitle('Bar jumps');
% subplot(3,1,1)
% plot(panels_y_raw)
% ylabel('Voltage (V)');xlabel('Time');
% subplot(3,1,2)
% plot(jumps);
% ylabel('Voltage difference (V)');xlabel('Time');
% 
% %plot the data from the yPanels and add lines of the previously determined
% %bar jumps
% subplot(3,1,3),
% plot(panels_y_raw)
% title('Bar jumps');
% xlabel('Time (frames)'); ylabel('Voltage (V)');
% hold on
% %add the bar jumps
% for i = 1:length(j)
%     plot([j(i) j(i)],[0 10],'r');
% end
% 

%% Fix the data after the jumps

%The x position of the panels and the FicTrac output are "ignorant" of the
%bar jumps. The x position of the bar will move with the angular position
%of the fly, but the coordinate system changes every time the bar jumps.
%We need to make sure the xpos and heading of the fly are corrected to take
%this coordinate change into account.

yVoltsBJ = panels_y_raw(j-1);
yVoltsAJ = panels_y_raw(j+1);
ydiff = yVoltsAJ-yVoltsBJ; %this is the offset that I need to adjust by

%Correct the data after every jump except for the last
panelsUW = panels;
ft_yawUW = ft_yaw;
for i = 1:size(j)-1
    panelsUW(j(i)+1:j(i+1)) = (panels(j(i)+1:j(i+1)))+sum(ydiff(1:i));
    ft_yawUW(j(i)+1:j(i+1)) = (ft_yaw(j(i)+1:j(i+1)))+sum(ydiff(1:i));
end

%Correct the data after the last jump
panelsUW(j(end)+1:end) = panels(j(end)+1:end)+sum(ydiff(1:end));
ft_yawUW(j(end)+1:end) = ft_yaw(j(end)+1:end)+sum(ydiff(1:end));

%I now have to wrap this data to get it to be between 0 and 10 V.
correct_panels = panelsUW;
correct_fr_yaw = ft_yawUW;

for i = 1:size(panels)
    
    if panelsUW(i) > 10
        correct_panels(i) = panelsUW(i)-10;
    else
        correct_panels(i) = panelsUW(i);
    end
    
    if ft_yawUW(i) > 10
        correct_ft_yaw(i) = ft_yawUW(i)-10;
    else
        correct_ft_yaw(i) = ft_yawUW(i);
    end
    
end

% %check the wrapping graphically
% figure('Position', [100 100 1200 900]),
% subplot(2,1,1),
% plot(panels_y_raw,'r')
% hold on
% plot(panels,'k')
% plot(panelsUW,'b')
% plot(correct_panels,'g')
% ylabel('Voltage (V)');
% title('x panels voltage signals before and after wrapping');
% ylim([0 max(panelsUW)]);
% legend({'yPanels','Uncorrected panels','Unwrapped xPanels', 'Wrapped xPanels'});
% 
% subplot(2,1,2),
% plot(panels_y_raw,'r')
% hold on
% plot(ft_yaw,'k')
% plot(ft_yawUW,'b')
% plot(correct_ft_yaw,'g')
% xlabel('Time (frame)');ylabel('Voltage (V)');
% title('angular position voltage signals before and after wrapping');
% ylim([0 max(ft_yawUW)]);
% legend({'yPanels','Uncorrected heading','Unwrapped angular position', 'Wrapped angular position'});


%% Recovering the jump magnitudes and make sure they are correct

jumps = [-120,120,-120,-120,120,120,-120,120,120,-120];

panels_y_pos = round ((panels_y_raw  * 96) /10);
jumpPos = panels_y_pos(j+1)-panels_y_pos(j-1);
degJumps = wrapTo180(jumpPos*(360/96));

panels_x_pos = round ((correct_panels  * 96) /10);
jumpMag = panels_x_pos(j+1)-panels_x_pos(j-1);
degMag = wrapTo180(jumpMag*(360/96));

jumpMag2 = correct_ft_yaw(j+1)-correct_ft_yaw(j-1);
radMag2 = jumpMag2.* 2 .* pi ./ 10; %convert from voltage to radians
degMag2 = wrapTo180(rad2deg(radMag2)); %convert from radians to degrees and wrap 180


% figure('Position', [100 100 1200 900]),
% 
% subplot(1,3,1), plot(degJumps,'ro')
% ylim([-180 180]);
% title('Jump magnitude, taken from yPanels');
% ylabel('deg');xlabel('Trial #');%xlim([1 jumps(end)]);
% % compare that to the jump function we have stored
% hold on
% plot(jumps,'b')
% legend({'Jumps from Y data','Jump function used'});
% 
% %Check if the jump magnitude appears ok in the x panel data
% subplot(1,3,2), plot(degMag,'ro')
% ylim([-180 180]);
% title('Jump magnitude, taken from xPanels');
% ylabel('deg');xlabel('Trial #');
% % compare that to the jump function we have stored
% hold on
% plot(jumps,'b')
% legend({'Jumps from X data','Jump function used'});
% 
% %check if the jump magnitude appears ok in the angular position data
% subplot(1,3,3), plot(degMag2,'ro')
% ylim([-180 180]);
% title('Jump magnitude, taken from angular position');
% ylabel('deg');xlabel('Trial #');
% % compare that to the jump function we have stored
% hold on
% plot(jumps,'b')
% legend({'Jumps from angular position','Jump function used'});


%% Process and save the correct data

data.ficTracIntx = ft_for;
data.ficTracInty = ft_side;
data.ficTracAngularPosition = correct_ft_yaw;
data.uncorrectedFicTracAngularPosition = ft_yaw;
[smoothed] = singleTrialVelocityAnalysis9mm(data,settings.sampRate);

[angle] = process_panel_360(correct_panels, num_frames_x);
flyPosRad = smoothed.angularPosition;
motor_pos = process_motor_360(motor);

[t] = process_time(trial_time);


end
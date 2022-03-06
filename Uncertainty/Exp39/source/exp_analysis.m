%code to analyze experiment 39 (cue jump experiment imaging EB-DANs)

clear all; close all;

%% Load data

[file,path] = uigetfile('Z:\Wilson Lab\Mel\Experiments\Uncertainty\Exp39\data\');
load(fullfile(path,file))

warning('off','all');

%% Make directory to save plots

%Move to the analysis folder
cd(path)
%List the contents
contents = dir();
%if there isn't a 'plots' folder already, create one
if (contains([contents.name],'plots') == 0)
   mkdir(path,'\plots'); 
end
%List the contents of the 'plots' folder
cd([path,'\plots\'])


%% Assign fly ID

folderNames = dir(path(1:53));
flyNames = struct();
for folder = 1:length(folderNames)
    if (contains(folderNames(folder).name,'75B10') & ~contains(folderNames(folder).name,'txt'))
        flyNames(folder).name = folderNames(folder).name;
    end
end
%Remove empty rows
flyNames = flyNames(~cellfun(@isempty,{flyNames.name}));

%Assign fly number
for fly = 1:length(flyNames)
    if strcmp(flyNames(fly).name,path(54:end-10))
        fly_ID = fly;
    end
end

%% Determine when the stimuli are on

panels_on = data.fr_y_ds>7;
wind_on = data.wind_valve>2;

figure,
subplot(2,1,1)
plot(panels_on)
ylim([-1 2]);
title('Panels on');

subplot(2,1,2)
plot(wind_on)
ylim([-1 2]);
title('Wind on');


%% Determine stimulus configuration

if mode(panels_on(1:100)) == 1
    configuration = 1; %bar first
else
    configuration = 2; %wind first
end

%% Look for change in stimuli

%Find the frames where the stimuli change
panels_change = abs(diff(panels_on));
panels_change_frames = find(panels_change>0.5);
wind_change = abs(diff(wind_on));
wind_change_frames = find(wind_change>0.5);

%Conversion factors
sec_to_frames = length(data.dff{1})/data.time(end);
frames_to_sec = data.time(end)/length(data.dff{1});


%% Find the jump frames

if configuration == 1
   %Find the bar jump frames
   coded_bar_jump_frames = [floor(1200*length(data.dff{1})/data.time(end)),floor(1800*length(data.dff{1})/data.time(end)),floor(2400*length(data.dff{1})/data.time(end)),floor(3000*length(data.dff{1})/data.time(end))];
   %Find the wind jump frames
   coded_wind_jump_frames = [floor(1500*length(data.dff{1})/data.time(end)),floor(2100*length(data.dff{1})/data.time(end)),floor(2700*length(data.dff{1})/data.time(end)),floor(3300*length(data.dff{1})/data.time(end))];   
else
   %Find the wind jump frames
   coded_wind_jump_frames = [floor(1200*length(data.dff{1})/data.time(end)),floor(1800*length(data.dff{1})/data.time(end)),floor(2400*length(data.dff{1})/data.time(end)),floor(3000*length(data.dff{1})/data.time(end))];
   %Find the bar jump frames
   coded_bar_jump_frames = [floor(1500*length(data.dff{1})/data.time(end)),floor(2100*length(data.dff{1})/data.time(end)),floor(2700*length(data.dff{1})/data.time(end)),floor(3300*length(data.dff{1})/data.time(end))];   
end

%Correct bar jump frames
abs_diff_bar_signal = abs(diff(unwrap(deg2rad(data.panel_angle))));
for jump = 1:length(coded_bar_jump_frames)
    frame_segment = [coded_bar_jump_frames(jump)-100:coded_bar_jump_frames(jump)+100];
    [~,I_bar_jump_frames(jump)] = max(abs_diff_bar_signal(frame_segment));
    real_bar_jump_frames(jump) = coded_bar_jump_frames(jump) + I_bar_jump_frames(jump) - 100;
end%%
if fly_ID == 1
    real_bar_jump_frames(1) = floor(1189.4*length(data.dff{1})/data.time(end));            
    real_bar_jump_frames(2) = floor(1789.4*length(data.dff{1})/data.time(end));        
    real_bar_jump_frames(3) = floor(2389.4*length(data.dff{1})/data.time(end));    
    real_bar_jump_frames(4) = floor(2989.5*length(data.dff{1})/data.time(end));
end
real_bar_jump_sec = real_bar_jump_frames*data.time(end)/length(data.dff{1});



%Correct wind jump frames
motor_pos = wrapTo180(rad2deg(data.motor_pos));
abs_diff_wind_signal = abs(diff(unwrap(deg2rad(motor_pos))));
for jump = 1:length(coded_wind_jump_frames)
    frame_segment = [coded_wind_jump_frames(jump)-100:coded_wind_jump_frames(jump)+100];
    [~,I_wind_jump_frames(jump)] = max(abs_diff_wind_signal(frame_segment));
    real_wind_jump_frames(jump) = coded_wind_jump_frames(jump) + I_wind_jump_frames(jump) - 100;
end

real_wind_jump_sec = real_wind_jump_frames*data.time(end)/length(data.dff{1});
%correct for the flies for which the method didn't work well
if fly_ID == 1
    real_wind_jump_frames(1) = floor(1489.58*length(data.dff{1})/data.time(end));
    real_wind_jump_frames(2) = floor(2089.58*length(data.dff{1})/data.time(end));            
    real_wind_jump_frames(3) = floor(2689.52*length(data.dff{1})/data.time(end));        
    real_wind_jump_frames(4) = floor(3289.65*length(data.dff{1})/data.time(end));    
end

%% Analyze full experiment

figure('Position',[100 100 1800 1000]),

%Plot the stimuli's position
subplot(3,1,1)
bar_pos = wrapTo180(data.panel_angle);
[x_out_bar,bar_to_plot] = removeWrappedLines(data.time(panels_on),bar_pos(panels_on));
plot(x_out_bar,bar_to_plot,'.','LineWidth',1.5)
hold on
%add motor position
[x_out_motor,motor_to_plot] = removeWrappedLines(data.time(wind_on),motor_pos(wind_on)');
plot(x_out_motor,motor_to_plot,'.','LineWidth',1.5)
%add jumps
xline(real_wind_jump_sec(1),'color',[0.4940 0.1840 0.5560],'linewidth',4)
xline(real_bar_jump_sec(1),'color',[0.4660 0.6740 0.1880],'linewidth',4)
legend('Wind jumps','Bar jumps');
for jump = 2:4
   xline(real_wind_jump_sec(jump),'color',[0.4940 0.1840 0.5560],'linewidth',4,'handlevisibility','off')
   xline(real_bar_jump_sec(jump),'color',[0.4660 0.6740 0.1880],'linewidth',4,'handlevisibility','off')
end
title('Stimuli position');
legend('Bar position','Wind position')
ylim([-180 180]);
xlim([0 x_out_bar(end)]);
set(gca,'xticklabel',{[]})

%Plot the EB-DAN activity
subplot(3,1,2)
plot(data.time,data.dff{1})
xlim([1 data.time(end)]);
title('EB-DAN activity');

%Plot the fly's movement
subplot(3,1,3)
plot(data.time,data.total_mvt_ds)
title('Total movement (deg/s)');
xlim([1 data.time(end)]);

saveas(gcf,[path,'\plots\full_experiment.png']);

%% Plot stim position around the jumps

all_jump_frames = sort([real_bar_jump_frames,real_wind_jump_frames]);

for jump = 1:length(all_jump_frames)
    
    time_around_jump = data.time(all_jump_frames(jump)-5*sec_to_frames:all_jump_frames(jump)+5*sec_to_frames);
    time_around_jump = time_around_jump - data.time(all_jump_frames(jump));
    
    figure('Position',[100 100 1000 800]),
    %Plot the stimuli
    subplot(3,1,1)
    plot(time_around_jump, wrapTo180(data.panel_angle(all_jump_frames(jump)-5*sec_to_frames:all_jump_frames(jump)+5*sec_to_frames)))
    hold on
    plot(time_around_jump, wrapTo180(rad2deg(data.motor_pos(all_jump_frames(jump)-5*sec_to_frames:all_jump_frames(jump)+5*sec_to_frames))))
    line([0 0],[-180 180],'color','r');    
    ylim([-180 180]);
    legend('Bar position','Wind position');
    ylabel('Position (deg)');
    
    %Plot the EB-DAN activity
    subplot(3,1,2)
    plot(time_around_jump, data.dff{1}(all_jump_frames(jump)-5*sec_to_frames:all_jump_frames(jump)+5*sec_to_frames))
    hold on
    line([0 0],[0 1],'color','r');
    ylabel('DF/F');
    
    %Plot the fly's rotational speed
    subplot(3,1,3)
    plot(time_around_jump, abs(data.vel_yaw_ds(all_jump_frames(jump)-5*sec_to_frames:all_jump_frames(jump)+5*sec_to_frames)))
    hold on
    line([0 0],[0 150],'color','r');
    ylabel('Fly rotational speed (deg/s)');
    xlabel('Time (sec)');
    
    suptitle(['Jump #',num2str(jump)]);
    
    saveas(gcf,[path,'plots\around_jump_',num2str(jump)])
    
end

%% Run the prediction analysis, using each of the ROI to define the DF/F

for ROI = 1:length(data.dff) 
    
    if ROI == 1
        roi_name = 'EB';
    elseif ROI == 2
        roi_name = 'left_bulb';
    else
        roi_name = 'right_bulb';
    end
        
    %% Cross correlation between EB-DAN activity and fly movement
    
    [c,lags] = xcorr(data.total_mvt_ds,data.dff{ROI});
    
    figure('Position',[100 100 1400 800]),
    stem(lags,c)
    hold on
    
    %Find the maximum cross-corr lag
    [M, index] = max(c);
    max_corr_diff = lags(index);
    
    line([max_corr_diff, max_corr_diff], [0 max(c)], 'color','r')
    title(['Dff lags movement by ',num2str(max_corr_diff),' frames']);
    
    saveas(gcf,[path,'plots\cross_corr_',roi_name,'.png']);
    
    %% Determine relationship between EB-DAN activity and movement in the three min preceding each jump
    
    %initialize coefficients vectors
    intercept = [];
    movement_coef = [];
    moving_true = [];
    adj_rs = [];
    
    for jump = 1:length(all_jump_frames)
        
        pre_jump_total_mvt = data.total_mvt_ds(all_jump_frames(jump)+max_corr_diff-180*sec_to_frames:all_jump_frames(jump)+max_corr_diff); %I'm looking 5 frames earlier into the movement
        pre_jump_mvt = abs(data.vel_yaw_ds(all_jump_frames(jump)+max_corr_diff-180*sec_to_frames:all_jump_frames(jump)+max_corr_diff)); %Im looking 5 frames earlier into the movement
        pre_jump_dff = data.dff{ROI}(all_jump_frames(jump)-180*sec_to_frames:all_jump_frames(jump));
        %create binary data for whether the animal is moving or not
        moving = pre_jump_total_mvt>20;
        %create data table
        all_data = table(pre_jump_mvt',pre_jump_dff',moving','VariableNames',{'movement','activity','moving'});
        all_data.moving = categorical(all_data.moving); %make the binary variable a categorical one
        %run model
        movement_model = fitlm(all_data,'activity ~ movement + moving');
        intercept = [intercept,movement_model.Coefficients.Estimate(1)];
        movement_coef = [movement_coef,movement_model.Coefficients.Estimate(2)];
        moving_true = [moving_true,movement_model.Coefficients.Estimate(3)];
        adj_rs = [adj_rs,movement_model.Rsquared.Adjusted];
        
    end
    
    %% Predict mean activity following jump based on movement
    
    %)1 ten second window
    
    for jump = 1:length(all_jump_frames)
        
        post_jump_total_mvt = data.total_mvt_ds(all_jump_frames(jump)+max_corr_diff:all_jump_frames(jump)+max_corr_diff+10*sec_to_frames);
        post_jump_mvt = abs(data.vel_yaw_ds(all_jump_frames(jump)+max_corr_diff:all_jump_frames(jump)+max_corr_diff+10*sec_to_frames));
        post_jump_moving = post_jump_total_mvt > 20;
        predicted_post_jump_dff = intercept(jump) + movement_coef(jump)*post_jump_mvt + moving_true(jump)*post_jump_moving;
        mean_predicted_post_jump_dff(jump) = mean(predicted_post_jump_dff);
        real_post_jump_dff = data.dff{ROI}(all_jump_frames(jump):all_jump_frames(jump)+10*sec_to_frames);
        mean_real_post_jump_dff(jump) = mean(real_post_jump_dff);
        
        %Plot to compare
        figure,
        plot(predicted_post_jump_dff)
        hold on
        plot(real_post_jump_dff)
        legend('predicted','real');
        ylabel('DF/F');
        title(['Rsquared = ',num2str(adj_rs(jump))]);
        saveas(gcf,[path,'plots\predicted_vs_real_jump',num2str(jump),'_activity_10_sec.png']);
        
    end
    
    %combine the data
    mean_post_jump_activity = [mean_predicted_post_jump_dff;mean_real_post_jump_dff];
    post_jump_activity = [predicted_post_jump_dff,real_post_jump_dff];
    
    %Plot
    figure('Position',[100 100 1000 800]),
    subplot(1,2,1)
    plot(mean_post_jump_activity,'-o','color',[.5 .5 .5])
    hold on
    plot(mean(mean_post_jump_activity,2),'-ko','linewidth',2)
    xlim([0 3]);
    xticks([1 2]);
    xticklabels({'predicted','real'});
    title('Post jump EB-DAN activity');
    ylabel('DF/F');
    
    subplot(1,2,2)
    group = [repelem(1,1,length(predicted_post_jump_dff)),repelem(2,1,length(real_post_jump_dff))];
    boxplot(post_jump_activity,group)
    hold on
    plot([1 2],[mean(predicted_post_jump_dff),mean(real_post_jump_dff)], 'dk')
    xlim([0 3]);
    xticks([1 2]);
    xticklabels({'predicted','real'});
    title('Post jump EB-DAN activity');
    ylabel('DF/F');
    
    suptitle('Ten seconds post jump');
    
    saveas(gcf,[path,'plots\predicted_vs_real_',roi_name,'_activity_10_sec.png']);
    
    %% Use a test period of 10 sec before the 3 min used for the model to compare to the 10 sec following the jump
  
    for jump = 1:length(all_jump_frames)
       
        test_total_mvt = data.total_mvt_ds(all_jump_frames(jump)+max_corr_diff-190*sec_to_frames:all_jump_frames(jump)+max_corr_diff-180*sec_to_frames);
        test_mvt = abs(data.vel_yaw_ds(all_jump_frames(jump)+max_corr_diff-190*sec_to_frames:all_jump_frames(jump)+max_corr_diff-180*sec_to_frames));
        test_moving = test_total_mvt > 20;
        predicted_test_dff = intercept(jump) + movement_coef(jump)*test_mvt + moving_true(jump)*test_moving;
        mean_predicted_test_dff(jump) = mean(predicted_test_dff);
        real_test_dff = data.dff{ROI}(all_jump_frames(jump)-190*sec_to_frames:all_jump_frames(jump)-180*sec_to_frames);
        mean_real_test_dff(jump) = mean(real_test_dff);
        
        %Plot to compare
        figure,
        plot(predicted_test_dff)
        hold on
        plot(real_test_dff)
        legend('predicted','real');
        ylabel('DF/F');
        ylim([-0.25 1]);
        saveas(gcf,[path,'plots\predicted_vs_real_jump',num2str(jump),'_test_activity_10_sec.png']);
        
    end
    
    %combine the data
    mean_test_activity = [mean_predicted_test_dff;mean_real_test_dff];
    test_activity = [predicted_test_dff,real_test_dff];
    
    %Plot
    figure('Position',[100 100 1000 800]),
    subplot(1,2,1)
    plot(mean_test_activity,'-o','color',[.5 .5 .5])
    hold on
    plot(mean(mean_test_activity,2),'-ko','linewidth',2)
    xlim([0 3]);
    xticks([1 2]);
    xticklabels({'predicted','real'});
    title('Post jump EB-DAN activity');
    ylabel('DF/F');
    
    subplot(1,2,2)
    group = [repelem(1,1,length(predicted_test_dff)),repelem(2,1,length(real_test_dff))];
    boxplot(test_activity,group)
    hold on
    plot([1 2],[mean(predicted_test_dff),mean(real_test_dff)], 'dk') %this seems to be wrong, as it isn't stored per jump? see for loop above
    xlim([0 3]);
    xticks([1 2]);
    ylim([-0.25 1]);
    xticklabels({'predicted','real'});
    title('Post jump EB-DAN activity');
    ylabel('DF/F');
    
    suptitle('Ten seconds test period');
    
    saveas(gcf,[path,'plots\predicted_vs_real_',roi_name,'_test_activity_10_sec.png']);
    
    %% Repeat using 2 sec as a window
    
    for jump = 1:length(all_jump_frames)
        
        post_jump_total_mvt = data.total_mvt_ds(all_jump_frames(jump)+max_corr_diff:all_jump_frames(jump)+max_corr_diff+2*sec_to_frames);
        post_jump_mvt = abs(data.vel_yaw_ds(all_jump_frames(jump)+max_corr_diff:all_jump_frames(jump)+max_corr_diff+2*sec_to_frames));
        post_jump_moving = post_jump_total_mvt > 20;
        predicted_post_jump_dff = intercept(jump) + movement_coef(jump)*post_jump_mvt + moving_true(jump)*post_jump_moving;
        mean_predicted_post_jump_dff(jump) = mean(predicted_post_jump_dff);
        real_post_jump_dff = data.dff{ROI}(all_jump_frames(jump):all_jump_frames(jump)+2*sec_to_frames);
        mean_real_post_jump_dff(jump) = mean(real_post_jump_dff);
        
        %Plot to compare
        figure,
        plot(predicted_post_jump_dff)
        hold on
        plot(real_post_jump_dff)
        legend('predicted','real');
        ylabel('DF/F');
        title(['Rsquared = ',num2str(adj_rs(jump))]);
        
    end
    
    %combine the data
    mean_post_jump_activity = [mean_predicted_post_jump_dff;mean_real_post_jump_dff];
    post_jump_activity = [predicted_post_jump_dff,real_post_jump_dff];
    
    %Plot
    figure('Position',[100 100 1000 800]),
    subplot(1,2,1)
    plot(mean_post_jump_activity,'-o','color',[.5 .5 .5])
    hold on
    plot(mean(mean_post_jump_activity,2),'-ko','linewidth',2)
    xlim([0 3]);
    xticks([1 2]);
    xticklabels({'predicted','real'});
    title('Post jump EB-DAN activity');
    ylabel('DF/F');
    
    subplot(1,2,2)
    group = [repelem(1,1,length(predicted_post_jump_dff)),repelem(2,1,length(real_post_jump_dff))];
    boxplot(post_jump_activity,group)
    hold on
    plot([1 2],[mean(predicted_post_jump_dff),mean(real_post_jump_dff)], 'dk')
    xlim([0 3]);
    xticks([1 2]);
    xticklabels({'predicted','real'});
    title('Post jump EB-DAN activity');
    ylabel('DF/F');
    
    suptitle('Two seconds post jump');
    
    saveas(gcf,[path,'plots\predicted_vs_real_',roi_name,'_activity_2_sec.png']);
    
    
    %% One second after the jump
    
    for jump = 1:length(all_jump_frames)
        
        post_jump_total_mvt = data.total_mvt_ds(all_jump_frames(jump)+max_corr_diff:all_jump_frames(jump)+max_corr_diff+1*sec_to_frames);
        post_jump_mvt = abs(data.vel_yaw_ds(all_jump_frames(jump)+max_corr_diff:all_jump_frames(jump)+max_corr_diff+1*sec_to_frames));
        post_jump_moving = post_jump_total_mvt > 20;
        predicted_post_jump_dff = intercept(jump) + movement_coef(jump)*post_jump_mvt + moving_true(jump)*post_jump_moving;
        mean_predicted_post_jump_dff(jump) = mean(predicted_post_jump_dff);
        real_post_jump_dff = data.dff{ROI}(all_jump_frames(jump):all_jump_frames(jump)+1*sec_to_frames);
        mean_real_post_jump_dff(jump) = mean(real_post_jump_dff);
        
        %Plot to compare
        figure,
        plot(predicted_post_jump_dff)
        hold on
        plot(real_post_jump_dff)
        legend('predicted','real');
        ylabel('DF/F');
        title(['Rsquared = ',num2str(adj_rs(jump))]);
        
    end
    
    %combine the data
    mean_post_jump_activity = [mean_predicted_post_jump_dff;mean_real_post_jump_dff];
    post_jump_activity = [predicted_post_jump_dff,real_post_jump_dff];
    
    %Plot
    figure('Position',[100 100 1000 800]),
    subplot(1,2,1)
    plot(mean_post_jump_activity,'-o','color',[.5 .5 .5])
    hold on
    plot(mean(mean_post_jump_activity,2),'-ko','linewidth',2)
    xlim([0 3]);
    xticks([1 2]);
    xticklabels({'predicted','real'});
    title('Post jump EB-DAN activity');
    ylabel('DF/F');
    
    subplot(1,2,2)
    group = [repelem(1,1,length(predicted_post_jump_dff)),repelem(2,1,length(real_post_jump_dff))];
    boxplot(post_jump_activity,group)
    hold on
    plot([1 2],[mean(predicted_post_jump_dff),mean(real_post_jump_dff)], 'dk')
    xlim([0 3]);
    xticks([1 2]);
    xticklabels({'predicted','real'});
    title('Post jump EB-DAN activity');
    ylabel('DF/F');
    
    suptitle('One second post jump');
    
    saveas(gcf,[path,'plots\predicted_vs_real_',roi_name,'_activity_1_sec.png']);
    
end

close all; clear all;
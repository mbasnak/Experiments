% dataRad_angularPosition = data.ficTracAngularPosition .* 2 .* pi ./ 10;
% 
% 
% %% 1) Using the unwrapped data to find jumps
% 
% unwrapped_angularPosition = unwrap(dataRad_angularPosition);
% 
% %find jumps
% figure
% plot(abs(diff(unwrapped_angularPosition)))
% jump_frames = find(abs(diff(unwrapped_angularPosition))>1);
% 
% %plot heading and add jump lines
% figure,
% plot(unwrapped_angularPosition)
% hold on
% for jump = 1:length(jump_frames)
%     xline(jump_frames(jump),'r')
% end
% 
% %recover the jump sizes in degrees
% jump_sizes_u = rad2deg(unwrapped_angularPosition(jump_frames+1) - unwrapped_angularPosition(jump_frames)); 
% figure,
% subplot(2,1,1)
% plot(jump_sizes_u,'-ro')
% title('Bar jumps recovered using unwrapped data');
% 
% %load the hd5f data to make sure the jump sizes match what python gave them
% hdf5_to_read = 'Z:\Wilson Lab\Mel\Experiments\Uncertainty\Exp30\data\20210304_60D05_7f_fly4\ball\hdf5_Stimulus_jump_Closed_Loop_20210304_165030_sid_1_tid_1000001.hdf5';
% 
% hd5f_jump_size = h5read(hdf5_to_read,'/jump_size');
% subplot(2,1,2)
% plot(hd5f_jump_size)
% 
% %% 2) Using the downsampled data to find jumps
% 
% resampled_angularPosition = resample(unwrapped_angularPosition,25,4000);
% 
% all_frames_time = linspace(1,900,length(unwrapped_angularPosition));
% ds_time = linspace(1,900,length(resampled_angularPosition));
% 
% %overlay regular and downsampled data
% for jump = 1:length(jump_frames)
%     figure,
%     plot(all_frames_time(jump_frames(jump)-3200:jump_frames(jump)+3200),unwrapped_angularPosition(jump_frames(jump)-3200:jump_frames(jump)+3200),'-o')
%     hold on
%     plot(ds_time(floor(jump_frames(jump)*25/4000-20):floor(jump_frames(jump)*25/4000+20)),resampled_angularPosition(floor(jump_frames(jump)*25/4000-20):floor(jump_frames(jump)*25/4000+20)),'-o')
%     title(['Jump #',num2str(jump)]);
%     
%     %saveas(gcf,['Z:\Wilson Lab\Mel\Experiments\Uncertainty\Exp30\data\extra_plots\jump',num2str(jump),'_js_RS.png']);
% end
% 
% 
% %% 3) Testing 'downsample' instead
% 
% downsampled_angularPosition = downsample(unwrapped_angularPosition,160);
% 
% %overlay regular and downsampled data
% for jump = 1:length(jump_frames)
%     figure,
%     plot(all_frames_time(jump_frames(jump)-3200:jump_frames(jump)+3200),unwrapped_angularPosition(jump_frames(jump)-3200:jump_frames(jump)+3200),'-o')
%     hold on
%     plot(ds_time(floor(jump_frames(jump)*25/4000-20):floor(jump_frames(jump)*25/4000+20)),downsampled_angularPosition(floor(jump_frames(jump)*25/4000-20):floor(jump_frames(jump)*25/4000+20)),'-o')
%     title(['Jump #',num2str(jump)]);
%     
%     %saveas(gcf,['Z:\Wilson Lab\Mel\Experiments\Uncertainty\Exp30\data\extra_plots\jump',num2str(jump),'_js_DS.png']);
% end
% 
% %% Find jumps with downsampled data
% 
% figure
% plot(abs(diff(downsampled_angularPosition)))
% jump_frames = find(abs(diff(downsampled_angularPosition))>1);
% 
% %plot heading and add jump lines
% figure,
% plot(downsampled_angularPosition)
% hold on
% for jump = 1:length(jump_frames)
%     xline(jump_frames(jump),'r')
% end
% 
% %recover the jump sizes in degrees
% jump_sizes_d = rad2deg(downsampled_angularPosition(jump_frames+1) - downsampled_angularPosition(jump_frames)); 
% figure,
% plot(jump_sizes_u,'-bo')
% hold on
% plot(jump_sizes_d,'-ro')
% title('Bar jumps recovered');
% legend('unwrapped','downsampled');
% 
% 
% %% Look at the jumps when smoothing data
% 
% smoothed.angularPosition = smoothdata(downsampled_angularPosition,'rlowess',25);
% 
% 
% %overlay regular and smoothed data
% for jump = 1:length(jump_frames)
%     figure,
%     plot(all_frames_time(jump_frames(jump)-3200:jump_frames(jump)+3200),unwrapped_angularPosition(jump_frames(jump)-3200:jump_frames(jump)+3200),'-o')
%     hold on
%     plot(ds_time(floor(jump_frames(jump)*25/4000-20):floor(jump_frames(jump)*25/4000+20)),smoothed.angularPosition(floor(jump_frames(jump)*25/4000-20):floor(jump_frames(jump)*25/4000+20)),'-o')
%     title(['Jump #',num2str(jump)]);
%     
%     %saveas(gcf,['Z:\Wilson Lab\Mel\Experiments\Uncertainty\Exp30\data\extra_plots\jump',num2str(jump),'_js_DS.png']);
% end
% 
% 
% %% Plot panel output
% 
% figure,
% subplot(7,1,1)
% plot(rawData)
% 
% subplot(7,1,2)
% plot(inv_data)
% 
% subplot(7,1,3)
% plot(deg_data)
% 
% subplot(7,1,4)
% plot(unwrapped_data)
% 
% subplot(7,1,5)
% plot(shifted_data)
% 
% subplot(7,1,6)
% plot(ds_data)
% 
% subplot(7,1,7)
% plot(stim_pos)
% 
% 
% %% Plot and find jump frames using this bar signal
% 
% figure,
% plot(diff(unwrap(deg2rad(stim_pos))))
% 
% jump_frames = find(abs(diff(unwrap(deg2rad(stim_pos))))>1)
% 
% jump_sizes = wrapTo180(stim_pos(jump_frames+1)-stim_pos(jump_frames));
% 
% figure,
% plot(jump_sizes,'-o')
% 
% 
% %% Plot and find jump frames using this bar signal ds
% 
% figure,
% plot(diff(unwrap(deg2rad(stim_pos_ds))))
% 
% jump_frames = find(abs(diff(unwrap(deg2rad(stim_pos_ds))))>1.5);
% 
% jump_sizes = wrapTo180(stim_pos_ds(jump_frames+1)-stim_pos_ds(jump_frames));
% 
% figure,
% plot(jump_sizes,'-o')
% 
% 
% %% Plot bar position around jumps
% 
% figure,
% plot(bar_position)
% for jump = 1:length(ds_jump_frames)
%     xline(ds_jump_frames(jump),'r')
% end
% 
% for jump = 1:length(ds_jump_frames)
%    figure,
%    plot(bar_position(ds_jump_frames(jump)-10:ds_jump_frames(jump)+10),'-o')
%    hold on
%    xline(11,'r')    
% end
% 
% jump_size = wrapTo180(bar_position(ds_jump_frames+1)-bar_position(ds_jump_frames));
% figure,plot(jump_size,'-ro')
% 
% 
% %% Find time it takes for offset to recover
% 
% figure
% plot(offset(1,:),'-o')
% hold on
% xline(26,'r')
% ylim([-180 180]);
% xlim([1 51]);

%% 
panels_rad = panels*2*pi/10;
unwrapped_panels = unwrap(panels_rad);
diff_panels = abs(diff(unwrapped_panels));

figure,
plot(diff_panels)

%remove last jump
jump_frames = find(diff_panels>1.5);
jump_frames(end) = [];
jump_frames(end) = [];

jump_sizes = rad2deg(circ_dist(panels_rad(jump_frames),panels_rad(jump_frames+1)));
figure,
plot(jump_sizes,'-ro');
ylabel('Bar jump size');
xlabel('Jump #');

%% downsampled data

ds_data = downsample(panels_rad, 4000/25);

unwrapped_panels2 = unwrap(ds_data);
diff_panels2 = abs(diff(unwrapped_panels2));

figure,
plot(diff_panels2)

%remove last jump
jump_frames2 = find(diff_panels2>1.5);

jump_sizes2 = rad2deg(circ_dist(ds_data(jump_frames2),ds_data(jump_frames2+1)));
figure,
plot(jump_sizes2,'-ro');
ylabel('Bar jump size');
xlabel('Jump #');

%% downsample again, to match the imaging rate

stim_pos_ds = ds_data(round(linspace(1, length(ds_data), volumes)));

unwrapped_panels3 = unwrap(stim_pos_ds);
diff_panels3 = abs(diff(unwrapped_panels3));

figure,
plot(diff_panels3)

%remove last jump
jump_frames3 = find(diff_panels3>1.5);

jump_sizes3 = rad2deg(circ_dist(stim_pos_ds(jump_frames3),stim_pos_ds(jump_frames3+1)));
figure,
plot(jump_sizes3,'-ro');
ylabel('Bar jump size');
xlabel('Jump #');

%% Overlay

figure,
plot(jump_sizes,'-bo');
hold on
plot(jump_sizes2,'-ko');
plot(jump_sizes3,'-go');
ylabel('Bar jump size');
xlabel('Jump #');
legend('Raw jump sizes','Ds jump sizes','Ds to imaging jump sizes');

%% look at individual jumps

%create time vectors
full_time = linspace(1,trial_time(end),length(panels_rad));
ds_time = linspace(1,trial_time(end),length(ds_data));
ds_ds_time = linspace(1,trial_time(end),length(stim_pos_ds));

for jump = [5,6,14]
    figure,
    plot(full_time(jump_frames(jump)-8000:jump_frames(jump)+8000),panels_rad(jump_frames(jump)-8000:jump_frames(jump)+8000),'b','linewidth',2)
    hold on
    %add the pre and post jump timepoints
    plot(full_time(jump_frames(jump)),panels_rad(jump_frames(jump)),'bo','linewidth',2)
    plot(full_time(jump_frames(jump)+1),panels_rad(jump_frames(jump)+1),'bo','linewidth',2)
    %add the jump size
    full_jump_size = rad2deg(circ_dist(panels_rad(jump_frames(jump)),panels_rad(jump_frames(jump)+1)));
    text(full_time(jump_frames(jump)+1000),pi*3/2,['Raw jump size = ',num2str(full_jump_size)]);
    
    %add the pre and post jump timepoints
    plot(ds_time(jump_frames2(jump)),ds_data(jump_frames2(jump)),'ko','linewidth',2)
    plot(ds_time(jump_frames2(jump)+1),ds_data(jump_frames2(jump)+1),'ko','linewidth',2)
    plot(ds_time(jump_frames2(jump)-50:jump_frames2(jump)+50),ds_data(jump_frames2(jump)-50:jump_frames2(jump)+50),'k','linewidth',2)
    %add the jump size
    ds_jump_size = rad2deg(circ_dist(ds_data(jump_frames2(jump)),ds_data(jump_frames2(jump)+1)));
    text(ds_time(jump_frames2(jump)+10),pi,['Ds jump size = ',num2str(ds_jump_size)],'color','b');
    
    %add the pre and post jump timepoints
    plot(ds_ds_time(jump_frames3(jump)),stim_pos_ds(jump_frames3(jump)),'go','linewidth',2)
    plot(ds_ds_time(jump_frames3(jump)+1),stim_pos_ds(jump_frames3(jump)+1),'go','linewidth',2)
    plot(ds_ds_time(jump_frames3(jump)-50:jump_frames3(jump)+50),stim_pos_ds(jump_frames3(jump)-50:jump_frames3(jump)+50),'g','linewidth',2)
    %add the jump size
    ds_ds_jump_size = rad2deg(circ_dist(stim_pos_ds(jump_frames3(jump)),stim_pos_ds(jump_frames3(jump)+1)));
    text(ds_ds_time(jump_frames3(jump)+2),pi*2/3,['Ds imaging jump size = ',num2str(ds_ds_jump_size)],'color','g');   
    
    
    %add plot aesthetics
    xline(ds_time(jump_frames2(jump)),'r','linewidth',2);
    ylim([0 2*pi]);
    xlim([ds_time(jump_frames2(jump)-50) ds_time(jump_frames2(jump)+50)]);
    title(['Jump #',num2str(jump)]);
    xlabel('Time (sec)');
    ylabel('Bar position (rad)');
end


%% Draw removing the jumps

%create derivative vectors
diff_data = abs(diff(panels_rad));
diff_ds_data = abs(diff(ds_data));
diff_ds_ds_data = abs(diff(stim_pos_ds));



for jump = [5,6,14]
    figure,
    %find jumps larger than 180 in data (i.e., where the data wraps around)
    plot_time = jump_frames(jump)-8000:jump_frames(jump)+8000;
    jump_data = find(diff_data(plot_time)>100*2*pi/360);
    if length(jump_data)>1
        plot(full_time(plot_time(1:jump_data(1))),panels_rad(plot_time(1:jump_data(1))),'-b','linewidth',2)
        hold on
        for wrap_point = 1:length(jump_data)-1
            plot(full_time(plot_time(jump_data(wrap_point)+1:jump_data(wrap_point+1))),panels_rad(plot_time(jump_data(wrap_point)+1:jump_data(wrap_point+1))),'-b','linewidth',2)
        end
        plot(full_time(plot_time(jump_data(wrap_point(end)+1):end)),panels_rad(plot_time(jump_data(wrap_point(end)+1):end)),'-b','linewidth',2)
    else
        plot(full_time(plot_time),panels_rad(plot_time),'-b','linewidth',2)
    end
%     hold on
%     %add the pre and post jump timepoints
%     plot(full_time(jump_frames(jump)),panels_rad(jump_frames(jump)),'bo','linewidth',2)
%     plot(full_time(jump_frames(jump)+1),panels_rad(jump_frames(jump)+1),'bo','linewidth',2)
%     %add the jump size
%     full_jump_size = rad2deg(circ_dist(panels_rad(jump_frames(jump)),panels_rad(jump_frames(jump)+1)));
%     text(full_time(jump_frames(jump)+1000),pi*3/2,['Raw jump size = ',num2str(full_jump_size)]);
    
%     %add the pre and post jump timepoints
%     plot(ds_time(jump_frames2(jump)),ds_data(jump_frames2(jump)),'ko','linewidth',2)
%     plot(ds_time(jump_frames2(jump)+1),ds_data(jump_frames2(jump)+1),'ko','linewidth',2)
%     plot(ds_time(jump_frames2(jump)-50:jump_frames2(jump)+50),ds_data(jump_frames2(jump)-50:jump_frames2(jump)+50),'k','linewidth',2)
%     %add the jump size
%     ds_jump_size = rad2deg(circ_dist(ds_data(jump_frames2(jump)),ds_data(jump_frames2(jump)+1)));
%     text(ds_time(jump_frames2(jump)+10),pi,['Ds jump size = ',num2str(ds_jump_size)],'color','b');
%     
%     %add the pre and post jump timepoints
%     plot(ds_ds_time(jump_frames3(jump)),stim_pos_ds(jump_frames3(jump)),'go','linewidth',2)
%     plot(ds_ds_time(jump_frames3(jump)+1),stim_pos_ds(jump_frames3(jump)+1),'go','linewidth',2)
%     plot(ds_ds_time(jump_frames3(jump)-50:jump_frames3(jump)+50),stim_pos_ds(jump_frames3(jump)-50:jump_frames3(jump)+50),'g','linewidth',2)
%     %add the jump size
%     ds_ds_jump_size = rad2deg(circ_dist(stim_pos_ds(jump_frames3(jump)),stim_pos_ds(jump_frames3(jump)+1)));
%     text(ds_ds_time(jump_frames3(jump)+2),pi*2/3,['Ds imaging jump size = ',num2str(ds_ds_jump_size)],'color','g');   
%     
    
    %add plot aesthetics
    xline(ds_time(jump_frames2(jump)),'r','linewidth',2);
    ylim([0 2*pi]);
    xlim([ds_time(jump_frames2(jump)-50) ds_time(jump_frames2(jump)+50)]);
    title(['Jump #',num2str(jump)]);
    xlabel('Time (sec)');
    ylabel('Bar position (rad)');
end
function run_gain_change_analysis()

%function to run the analysis for the closed-loop bouts of the change in
%contrast experiment in all the flies of the dataset



%Get the path for each fly
parentDir = 'Z:\Wilson Lab\Mel\Experiments\Uncertainty\Exp28\data';
folderNames = dir(parentDir);

for content = 1:length(folderNames)
    if contains(folderNames(content).name,'60D05')
        flyData{content} = [folderNames(content).folder,'\',folderNames(content).name];
    end
end

%remove empty cells
data_dirs = flyData(~cellfun(@isempty,flyData));

%determine the session id for this trial type from the session_info.mat
%file, load the data and plot
for fly = 1:length(data_dirs)
    
    clearvars -except data_dirs fly
    load([data_dirs{fly},'\sessions_info.mat']);
    sid = sessions_info.gain_change;
    
    %get the contents of the fly folder
    fly_files = dir([data_dirs{fly},'\analysis']);
    
    %determine which content belongs to the sid we extracted
    for file = 1:length(fly_files)
        if contains(fly_files(file).name,['sid_',num2str(sid),'_'])
            %load the data
            load([fly_files(file).folder,'\',fly_files(file).name])
            
            path = [fly_files(file).folder,'\'];
            
            %% Determine changes in gain
            
            %determine initial gain
            %get hdf5 files in directory
            pause(1);
            file_names = dir(fullfile([path(1:end-9),'ball\'],'*hdf5'));
            for file = 1:length(file_names)
                if contains(file_names(file).name,['sid_',num2str(sid)])
                    hdf5_file_to_read = fullfile(file_names(file).folder,file_names(file).name);
                end
            end
            
            gain_yaw = double(h5read(hdf5_file_to_read,'/gain_yaw'));
            %downsample to match data length
            gain_yaw_ds = resample(gain_yaw,length(data.time),length(gain_yaw));
            gain_yaw_ds(gain_yaw_ds<0) = -1;
            gain_yaw_ds(gain_yaw_ds>0) = 1;
            
            %determine gain changes
            gain_changes = find(abs(diff(gain_yaw_ds))>0.5);
            gain_changes = gain_changes(1:2);
            
            %% Set block limits
            
            blockLimits{1} = [1,gain_changes(1)-1];
            blockLimits{2} = [gain_changes(1),gain_changes(2)];
            blockLimits{3} = [gain_changes(2)+1,length(data.time)];
            
            %determine gain per block
            gain_per_block = [1,-1,1];
            
            
            %set color palette based on gain
            for block = 1:3
                if gain_per_block(block) == 1
                    color_gradient{block} = [.4 0.1 0.1];
                else
                    color_gradient{block} = [.1 0.4 0.2];
                end
            end
            
            %% Plot heatmap and offset
            
            %Plot
            % Plot the heatmap of EPG activity
            figure('Position',[100 100 1400 800]),
            subplot(4,5,[1 5])
            imagesc(data.dff_matrix)
            colormap(flipud(gray))
            %colormap(gray)
            hold on
            %add the changes in stim
            for change = 1:length(gain_changes)
                line([gain_changes(change) gain_changes(change)], [0 17], 'LineWidth', 3, 'color', [0, 0.5, 0]);
            end
            yticks(1:2:16);
            yticklabels({'8L','6L','4L','2L','1R','3R','5R','7R'});
            ylabel('PB glomerulus');
            title('EPG activity in the PB');
            set(gca,'XTickLabel',[]);
            legend('Change in stimulus');
            
            % Plot the bar position and the EPG phase
            subplot(4,5,[6 10])
            %Get bar position to plot
            bar_position = wrapTo180(data.panel_angle);
            % change_bar_pos = abs([0;diff(smooth(bar_position))]);
            % bar_pos_to_plot = smooth(bar_position);
            % bar_pos_to_plot(change_bar_pos>40==1) = NaN;
            [x_out_bar, bar_pos_to_plot] = removeWrappedLines(data.time,bar_position);
            plot(x_out_bar,bar_pos_to_plot,'LineWidth',1.5,'color',[0.2 0.6 0.7])
            hold on
            %Get EPG phase to plot
            %I'm now going to negate the phase, since I'm plotting heading instead of
            %bar position, to the bump moves in the other direction
            phase = wrapTo180(rad2deg(data.dff_pva));
            % change_phase = abs([0;diff(smooth(phase))]);
            % phase_to_plot = smooth(phase);
            % phase_to_plot(change_phase>40==1) = NaN;
            [x_out_phase, phase_to_plot] = removeWrappedLines(data.time,phase);
            plot(x_out_phase,phase_to_plot,'LineWidth',1.5,'color',[0.9 0.3 0.4])
            %add the changes in stim
            for change = 1:length(gain_changes)
                line([data.time(gain_changes(change)) data.time(gain_changes(change))], [-180 180], 'LineWidth', 2, 'color',[0, 0.5, 0]);
            end
            legend('Bar position', 'EPG phase');
            title('EPG phase and bar position');
            ylim([-180, 180]);
            xlim([0,data.time(end)]);
            ylabel('Deg');
            set(gca,'XTickLabel',[]);
            
            % Plot the offset
            subplot(4,5,[11 15])
            %Get offset to plot
            offset = wrapTo180(data.bar_offset);
            % change_offset = abs([0;diff(smooth(offset))]);
            % offset_to_plot = smooth(offset);
            % offset_to_plot(change_offset>30==1) = NaN;
            [x_out_offset, offset_to_plot] = removeWrappedLines(data.time,offset);
            plot(x_out_offset,offset_to_plot,'LineWidth',1.5,'color','k')
            %add the changes in stim
            for change = 1:length(gain_changes)
                line([data.time(gain_changes(change)) data.time(gain_changes(change))], [-180 180], 'LineWidth', 2, 'color', [0, 0.5, 0]);
            end
            ylim([-180 180]);
            title('Bar offset');
            ylabel('Deg'); xlabel('Time (sec)');
            
            % Polar histograms of offset
            % Color histograms acording to the intensity level of the bar, using the
            % vector Intensities
            subplot(4,5,16)
            polarhistogram(deg2rad(offset(blockLimits{1}(1):blockLimits{1}(2))),15,'FaceColor',color_gradient{1})
            set(gca,'ThetaZeroLocation','top',...
                'ThetaDir','counterclockwise');
            title({'Offset distribution';'normal gain'});
            rticklabels({})
            
            subplot(4,5,18)
            polarhistogram(deg2rad(offset(blockLimits{2}(1):blockLimits{2}(2))),15,'FaceColor',color_gradient{2})
            set(gca,'ThetaZeroLocation','top',...
                'ThetaDir','counterclockwise');
            title({'Offset distribution';'inverted gain'});
            rticklabels({})
            
            subplot(4,5,20)
            polarhistogram(deg2rad(offset(blockLimits{3}(1):blockLimits{3}(2))),15,'FaceColor',color_gradient{3})
            set(gca,'ThetaZeroLocation','top',...
                'ThetaDir','counterclockwise');
            title({'Offset distribution';'normal gain'});
            rticklabels({})
            
            %save figure
            saveas(gcf,[path,'plots\gainChangeHeatmapAndBarOffset.png']);
            
            %% Heatmap using heading instead of bar position
            
            %Plot
            % Plot the heatmap of EPG activity
            figure('Position',[100 100 1400 800]),
            subplot(4,5,[1 5])
            %I will now flip the dff matrix, since I'm plotting with respect to heading
            imagesc(flip(data.dff_matrix))
            colormap(flipud(gray))
            %colormap(gray)
            hold on
            %add the changes in stim
            for change = 1:length(gain_changes)
                line([gain_changes(change) gain_changes(change)], [0 17], 'LineWidth', 3, 'color', [0, 0.5, 0]);
            end
            yticks(1:2:16);
            yticklabels({'8R','6R','4R','2R','1L','3L','5L','7L'});
            ylabel('PB glomerulus');
            title('EPG activity in the PB');
            set(gca,'XTickLabel',[]);
            legend('Change in stimulus');
            
            % Plot the bar position and the EPG phase
            subplot(4,5,[6 10])
            %Get bar position to plot
            heading = wrapTo180(data.heading_deg);
            [x_out_heading, heading_to_plot] = removeWrappedLines(data.time,heading);
            plot(x_out_heading,heading_to_plot,'LineWidth',1.5,'color',[0.6 0.3 0.8])
            hold on
            %Get EPG phase to plot
            %I'm now going to negate the phase, since I'm plotting heading instead of
            %bar position, to the bump moves in the other direction
            phase = wrapTo180(rad2deg(-data.dff_pva));
            [x_out_phase, phase_to_plot] = removeWrappedLines(data.time,phase);
            plot(x_out_phase,phase_to_plot,'LineWidth',1.5,'color',[0.9 0.3 0.4])
            %add the changes in stim
            for change = 1:length(gain_changes)
                line([data.time(gain_changes(change)) data.time(gain_changes(change))], [-180 180], 'LineWidth', 2, 'color',[0, 0.5, 0]);
            end
            legend('Fly heading', 'EPG phase');
            title('EPG phase and bar position');
            ylim([-180, 180]);
            xlim([0,data.time(end)]);
            ylabel('Deg');
            set(gca,'XTickLabel',[]);
            
            % Plot the offset
            subplot(4,5,[11 15])
            %Get offset to plot
            heading_offset = wrapTo180(data.heading_offset);
            [x_out_heading_offset, heading_offset_to_plot] = removeWrappedLines(data.time,heading_offset);
            plot(x_out_heading_offset,heading_offset_to_plot,'LineWidth',1.5,'color','k')
            %add the changes in stim
            for change = 1:length(gain_changes)
                line([data.time(gain_changes(change)) data.time(gain_changes(change))], [-180 180], 'LineWidth', 2, 'color', [0, 0.5, 0]);
            end
            ylim([-180 180]);
            title('Heading offset');
            ylabel('Deg'); xlabel('Time (sec)');
            
            % Polar histograms of offset
            % Color histograms acording to the intensity level of the bar, using the
            % vector Intensities
            subplot(4,5,16)
            polarhistogram(deg2rad(heading_offset(blockLimits{1}(1):blockLimits{1}(2))),15,'FaceColor',color_gradient{1})
            set(gca,'ThetaZeroLocation','top',...
                'ThetaDir','counterclockwise');
            title({'Offset distribution';'normal gain'});
            rticklabels({})
            
            subplot(4,5,18)
            polarhistogram(deg2rad(heading_offset(blockLimits{2}(1):blockLimits{2}(2))),15,'FaceColor',color_gradient{2})
            set(gca,'ThetaZeroLocation','top',...
                'ThetaDir','counterclockwise');
            title({'Offset distribution';'inverted gain'});
            rticklabels({})
            
            subplot(4,5,20)
            polarhistogram(deg2rad(heading_offset(blockLimits{3}(1):blockLimits{3}(2))),15,'FaceColor',color_gradient{3})
            set(gca,'ThetaZeroLocation','top',...
                'ThetaDir','counterclockwise');
            title({'Offset distribution';'normal gain'});
            rticklabels({})
            
            %save figure
            saveas(gcf,[path,'plots\gainChangeHeatmapAndHeadingOffset.png']);
            
            
            %% Offset stability
            
            %I will analyze offset stability by computing a rolling window of circular
            %standard deviation of offset and taking the inverse
            fcn = @(x) circ_std(x);
            heading_offset_variability = matlab.tall.movingWindow(fcn,50,deg2rad(heading_offset));
            smoothed_heading_offset_variability = smooth(heading_offset_variability,150,'rloess');
            
            figure('Position',[100 100 1600 400]),
            subplot(3,1,1)
            plot(x_out_heading_offset,heading_offset_to_plot,'k')
            hold on
            for change = 1:length(gain_changes)
                line([data.time(gain_changes(change)) data.time(gain_changes(change))], [-180 180], 'LineWidth', 2, 'color', [0, 0.5, 0]);
            end
            ylim([-180 180]);
            title('Heading offset');
            set(gca,'xticklabel',{[]});
            
            subplot(3,1,2)
            plot(data.time,heading_offset_variability,'k')
            hold on
            for change = 1:length(gain_changes)
                line([data.time(gain_changes(change)) data.time(gain_changes(change))], [0 2], 'LineWidth', 2, 'color', [0, 0.5, 0]);
            end
            title('Heading offset variability');
            set(gca,'xticklabel',{[]});
            
            subplot(3,1,3)
            plot(data.time,smoothed_heading_offset_variability,'k')
            hold on
            for change = 1:length(gain_changes)
                line([data.time(gain_changes(change)) data.time(gain_changes(change))], [0 2], 'LineWidth', 2, 'color', [0, 0.5, 0]);
            end
            title('Smoothed heading offset variability');
            xlabel('Time (sec)');
            
            %save figure
            saveas(gcf,[path,'plots\HeadingOffsetVariability.png']);
            
            %% Repeat for bar offset
            
            offset_variability = matlab.tall.movingWindow(fcn,50,deg2rad(offset));
            smoothed_bar_offset_variability = smooth(offset_variability,150,'rloess');
            
            
            figure('Position',[100 100 1600 400]),
            subplot(3,1,1)
            plot(x_out_offset,offset_to_plot,'k')
            hold on
            for change = 1:length(gain_changes)
                line([data.time(gain_changes(change)) data.time(gain_changes(change))], [-180 180], 'LineWidth', 2, 'color', [0, 0.5, 0]);
            end
            ylim([-180 180]);
            title('Bar offset');
            set(gca,'xticklabel',{[]});
            
            subplot(3,1,2)
            plot(data.time,offset_variability,'k')
            hold on
            for change = 1:length(gain_changes)
                line([data.time(gain_changes(change)) data.time(gain_changes(change))], [0 2], 'LineWidth', 2, 'color', [0, 0.5, 0]);
            end
            title('Bar offset variability');
            set(gca,'xticklabel',{[]});
            
            subplot(3,1,3)
            plot(data.time,smoothed_bar_offset_variability,'k')
            hold on
            for change = 1:length(gain_changes)
                line([data.time(gain_changes(change)) data.time(gain_changes(change))], [0 2], 'LineWidth', 2, 'color', [0, 0.5, 0]);
            end
            title('Smoothed bar offset variability');
            xlabel('Time (sec)');
            
            %save figure
            saveas(gcf,[path,'plots\BarOffsetVariability.png']);
            
            %% Determine the type of fly based on the ratio between bar and offset variability
            
            ratio_means = mean(offset_variability)/mean(heading_offset_variability);
            ratio_medians = median(offset_variability)/median(heading_offset_variability);
            
            all_ratios = [ratio_means,ratio_medians];
            
            figure,
            plot(all_ratios,'-ko','linewidth',2)
            hold on
            yline(1,'r','linewidth',2)
            text(0.2,1.2,'Type 1');
            text(0.2,0.8,'Type 2');
            ylim([0 2]);
            xlim([0 3]);
            xticks([1 2])
            xticklabels({'ratio of means','ratio of medians'})
            ylabel('Bar offset variability / heading offset variability');
            %save figure
            saveas(gcf,[path,'plots\FlyClassification.png']);
            
            %classify fly
            if mean(all_ratios) > 1
                type_of_fly = 1; %fly that learns the new mapping (and where heading offset is more stable)
            else
                type_of_fly = 2; %fly that ignores proprioceptive cues (and where bar offset is more stable)
            end
            
            
            %% Calculate and plot bump magnitude in time
            
            %compute bump magnitude as max-min
            allBumpMag = [];
            allZBumpMag = [];
            for block = 1:length(blockLimits)
                bump_mag{block} = data.bump_magnitude(:,blockLimits{block}(1):blockLimits{block}(2));
                allBumpMag = [allBumpMag,bump_mag{block}];
                %add the bump mag computed with the zscored epg data
                z_bump_mag{block} = max(data.mean_z_EB(:,blockLimits{block}(1):blockLimits{block}(2)))-min(data.mean_z_EB(:,blockLimits{block}(1):blockLimits{block}(2)));
                allZBumpMag = [allZBumpMag,z_bump_mag{block}];
            end
            
            figure('Position',[200 200 1600 600]),
            %plot EPG activity
            subplot(3,1,1)
            imagesc(flip(data.dff_matrix))
            colormap(flipud(gray))
            hold on
            %add the changes in stim
            for change = 1:length(gain_changes)
                line([gain_changes(change) gain_changes(change)], [0 17], 'LineWidth', 2, 'color', [0, 0.5, 0]);
            end
            yticks(1:2:16);
            yticklabels({'8L','6L','4L','2L','1R','3R','5R','7R'});
            ylabel('PB glomerulus');
            title('EPG activity in the PB');
            set(gca,'XTickLabel',[]);
            legend('Change in stimulus');
            
            %plot bump magnitude
            subplot(3,1,2)
            for block = 1:length(blockLimits)
                plot(data.time(blockLimits{block}(1):blockLimits{block}(2)),smooth(bump_mag{block}),'color',color_gradient{block})
                hold on
            end
            title('Bump magnitude');
            ylabel('Bump magnitude (max-min)');
            set(gca,'xticklabel',{[]});
            xlim([0 data.time(end)]);
            
            
            % Bump width at half max
            for timepoint = 1:size(data.mean_dff_EB,2)
                
                %linearly interpolate to have 1000 datapoints instead of 8
                interp_ex_data = interp1([1:8],data.mean_dff_EB(:,timepoint),[1:7/1000:8]);
                %Find the half max point
                half_max = (max(data.mean_dff_EB(:,timepoint))-min(data.mean_dff_EB(:,timepoint)))/2 + min(data.mean_dff_EB(:,timepoint));
                [ex_bump_mag_interp I_interp] = max(interp_ex_data);
                %Find in each half the index closest to the half max
                diff_data = abs(interp_ex_data-half_max);
                [sortedVals,indexes] = sort(diff_data);
                %remove all the indexes that are less than 175 datapoints away (i.e., a little more than 1 glomerulus) from
                %index(1)
                diff_indexes = abs(indexes-indexes(1));
                indexes(diff_indexes<175 & diff_indexes>0)=NaN;
                indexes = indexes(~isnan(indexes));
                two_indexes = [indexes(1), indexes(2)];
                I1 = min(two_indexes);
                I2 = max(two_indexes);
                if (all(two_indexes>I_interp) | all(two_indexes<I_interp))
                    half_max_w = I1+1000-I2;
                else
                    half_max_w = I2-I1;
                end
                %convert to EB tiles
                half_max_width(timepoint) = half_max_w*8/1001;
                
            end
            
            % Analyze in time
            allHalfWidth = [];
            for block = 1:length(blockLimits)
                width_half_max{block} = half_max_width(blockLimits{block}(1):blockLimits{block}(2));
                allHalfWidth = [allHalfWidth,width_half_max{block}];
            end
            
            %plot
            subplot(3,1,3)
            for block = 1:length(blockLimits)
                plot(data.time(blockLimits{block}(1):blockLimits{block}(2)),smooth(width_half_max{block}),'color',color_gradient{block})
                hold on
            end
            title('Bump half width');
            ylabel('Bump half width');
            xlabel('Time (sec)');
            xlim([0 data.time(end)]);
            
            %save figure
            saveas(gcf,[path,'plots\gainChangeBMandHWinTime.png']);
            
            %% Calculate bump width at half max with the zscored data
            
            for timepoint = 1:size(data.mean_z_EB,2)
                
                %linearly interpolate to have 1000 datapoints instead of 8
                interp_ex_data = interp1([1:8],data.mean_z_EB(:,timepoint),[1:7/1000:8]);
                %Find the half max point
                half_max = (max(data.mean_z_EB(:,timepoint))-min(data.mean_z_EB(:,timepoint)))/2 + min(data.mean_z_EB(:,timepoint));
                [ex_bump_mag_interp I_interp] = max(interp_ex_data);
                %Find in each half the index closest to the half max
                diff_data = abs(interp_ex_data-half_max);
                [sortedVals,indexes] = sort(diff_data);
                %remove all the indexes that are less than 175 datapoints away (i.e., a little more than 1 glomerulus) from
                %index(1)
                diff_indexes = abs(indexes-indexes(1));
                indexes(diff_indexes<175 & diff_indexes>0)=NaN;
                indexes = indexes(~isnan(indexes));
                two_indexes = [indexes(1), indexes(2)];
                I1 = min(two_indexes);
                I2 = max(two_indexes);
                if (all(two_indexes>I_interp) | all(two_indexes<I_interp))
                    half_max_w = I1+1000-I2;
                else
                    half_max_w = I2-I1;
                end
                %convert to EB tiles
                z_half_max_width(timepoint) = half_max_w*8/1001;
                
            end
            
            % Analyze in time
            
            allZHalfWidth = [];
            for block = 1:length(blockLimits)
                z_width_half_max{block} = z_half_max_width(blockLimits{block}(1):blockLimits{block}(2));
                allZHalfWidth = [allZHalfWidth,z_width_half_max{block}];
            end
            
            %% New heatmap plot including bump magnitude
            
            if type_of_fly == 1
                
                % Plot the heatmap of EPG activity
                figure('Position',[100 100 1400 800]),
                subplot(5,1,1)
                %I'm flipping the dff matrix for it to make sense along with the fly's
                %heading
                imagesc(flip(data.dff_matrix))
                colormap(flipud(gray))
                hold on
                %add the changes in stim
                for change = 1:length(gain_changes)
                    line([gain_changes(change) gain_changes(change)], [0 17], 'LineWidth', 2, 'color', [0, 0.5, 0]);
                end
                yticks(1:2:16);
                yticklabels({'8L','6L','4L','2L','1R','3R','5R','7R'});
                ylabel('PB glomerulus');
                title('EPG activity in the PB');
                set(gca,'XTickLabel',[]);
                legend('Change in stimulus');
                
                % Plot the bar position and the EPG phase
                subplot(5,1,2)
                %Get bar position to plot
                plot(x_out_heading,heading_to_plot,'color',[0.6 0.3 0.8],'LineWidth',1.5)
                hold on
                plot(x_out_phase,phase_to_plot,'color',[0.9 0.3 0.4],'LineWidth',1.5)
                %add the changes in stim
                for change = 1:length(gain_changes)
                    line([data.time(gain_changes(change)) data.time(gain_changes(change))], [-180 180], 'LineWidth', 2, 'color',  [0, 0.5, 0]);
                end
                legend('Fly heading', 'EPG phase');
                title('EPG phase and fly heading');
                ylim([-180, 180]);
                xlim([0,data.time(end)]);
                ylabel('Deg');
                set(gca,'XTickLabel',[]);
                
                % Plot the offset
                subplot(5,1,3)
                %Get offset to plot
                plot(x_out_heading_offset,heading_offset_to_plot,'LineWidth',1.5,'color','k')
                %add the changes in stim
                for change = 1:length(gain_changes)
                    line([data.time(gain_changes(change)) data.time(gain_changes(change))], [-180 180], 'LineWidth', 2, 'color',  [0, 0.5, 0]);
                end
                ylim([-180 180]);
                ylabel('Deg');
                set(gca,'XTickLabel',[]);
                title('Heading offset');
                
                %Plot the bump magnitude
                subplot(5,1,4)
                for block = 1:length(blockLimits)
                    plot(data.time(blockLimits{block}(1):blockLimits{block}(2)),smooth(bump_mag{block}),'color',color_gradient{block},'linewidth',1.5)
                    hold on
                end
                title('Bump magnitude');
                ylabel({'Bump magnitude';'(amplitude of Fourier component)'});
                set(gca,'xticklabel',{[]});
                xlim([0 data.time(end)]);
                
                subplot(5,1,5)
                for block = 1:length(blockLimits)
                    plot(data.time(blockLimits{block}(1):blockLimits{block}(2)),smooth(width_half_max{block}),'color',color_gradient{block},'linewidth',1.5)
                    hold on
                end
                title('Bump half width');
                ylabel('Bump half width');
                xlabel('Time (sec)');
                xlim([0 data.time(end)]);
                
                %save figure
                saveas(gcf,[path,'plots\HeadingOffsetWithBumpParameters.png']);
                
                
            else
                
                % Plot the heatmap of EPG activity
                figure('Position',[100 100 1400 800]),
                subplot(5,1,1)
                imagesc(data.dff_matrix)
                colormap(flipud(gray))
                %colormap(gray)
                hold on
                %add the changes in stim
                for change = 1:length(gain_changes)
                    line([gain_changes(change) gain_changes(change)], [0 17], 'LineWidth', 2, 'color', [0,0.5,0]);
                end
                yticks(1:2:16);
                yticklabels({'8L','6L','4L','2L','1R','3R','5R','7R'});
                ylabel('PB glomerulus');
                title('EPG activity in the PB');
                set(gca,'XTickLabel',[]);
                legend('Change in stimulus');
                
                % Plot the bar position and the EPG phase
                subplot(5,1,2)
                %Get bar position to plot
                bar_position = wrapTo180(data.panel_angle);
                plot(x_out_bar,bar_pos_to_plot,'color',[0.2 0.6 0.7],'LineWidth',1.5)
                hold on
                %Get EPG phase to plot
                %I'm now going to negate the phase, since I'm plotting heading instead of
                %bar position, to the bump moves in the other direction
                phase = wrapTo180(rad2deg(data.dff_pva));
                [x_out_phase,phase_to_plot] = removeWrappedLines(data.time,phase);
                plot(x_out_phase,phase_to_plot,'color',[0.9 0.3 0.4],'LineWidth',1.5)
                %add the changes in stim
                for change = 1:length(gain_changes)
                    line([data.time(gain_changes(change)) data.time(gain_changes(change))], [-180 180], 'LineWidth', 2, 'color', [0,0.5,0]);
                end
                legend('Bar position', 'EPG phase');
                title('EPG phase and bar position');
                ylim([-180, 180]);
                xlim([0,data.time(end)]);
                ylabel('Deg');
                set(gca,'XTickLabel',[]);
                
                % Plot the offset
                subplot(5,1,3)
                %Get offset to plot
                plot(x_out_offset,offset_to_plot,'LineWidth',1.5,'color','k')
                %add the changes in stim
                for change = 1:length(gain_changes)
                    line([data.time(gain_changes(change)) data.time(gain_changes(change))], [-180 180], 'LineWidth', 2, 'color', [0,0.5,0]);
                end
                ylim([-180 180]);
                title('Bar offset');
                ylabel('Deg'); xlabel('Time (sec)');
                
                %Plot the bump magnitude
                subplot(5,1,4)
                for block = 1:length(blockLimits)
                    plot(data.time(blockLimits{block}(1):blockLimits{block}(2)),smooth(bump_mag{block}),'color',color_gradient{block},'linewidth',1.5)
                    hold on
                end
                title('Bump magnitude');
                ylabel({'Bump magnitude';'(amplitude of Fourier component)'});
                set(gca,'xticklabel',{[]});
                xlim([0 data.time(end)]);
                
                subplot(5,1,5)
                for block = 1:length(blockLimits)
                    plot(data.time(blockLimits{block}(1):blockLimits{block}(2)),smooth(width_half_max{block}),'color',color_gradient{block},'linewidth',1.5)
                    hold on
                end
                title('Bump half width');
                ylabel('Bump half width');
                xlabel('Time (sec)');
                xlim([0 data.time(end)]);
                
                %save figure
                saveas(gcf,[path,'plots\BarOffsetWithBumpParameters.png']);
                
            end
            
            
            
            %%  Correlate BM and HW with heading or bar offset variability depending on the type of fly
            
            if type_of_fly == 1
                
                figure('Position',[100 100 1000 800]),
                
                subplot(2,2,1)
                binscatter(heading_offset_variability,smooth(allBumpMag))
                correlations = corrcoef(heading_offset_variability,smooth(allBumpMag));
                ylabel('Bump magnitude (max-min)');
                xlabel('Heading offset variability');
                text(1,2.5,['Corr ',num2str(round(correlations(2,1),2))]);
                title('Non-smoothed data');
                
                subplot(2,2,2)
                binscatter(heading_offset_variability,smooth(allHalfWidth))
                correlations = corrcoef(heading_offset_variability,smooth(allHalfWidth));
                ylabel('Bump half width');
                xlabel('Heading offset variability');
                text(1,6,['Corr ',num2str(round(correlations(2,1),2))]);
                title('Non-smoothed data');
                
                %with smoothed offset variability
                subplot(2,2,3)
                binscatter(smoothed_heading_offset_variability,smooth(allBumpMag))
                correlations = corrcoef(smoothed_heading_offset_variability,smooth(allBumpMag));
                ylabel('Bump magnitude (max-min)');
                xlabel('Heading offset variability');
                text(1,2.5,['Corr ',num2str(round(correlations(2,1),2))]);
                title('Smoothed data');
                
                subplot(2,2,4)
                binscatter(smoothed_heading_offset_variability,smooth(allHalfWidth))
                correlations = corrcoef(smoothed_heading_offset_variability,smooth(allHalfWidth));
                ylabel('Bump half width');
                xlabel('Heading offset variability');
                text(1,6,['Corr ',num2str(round(correlations(2,1),2))]);
                title('Smoothed data');
                
                
                
                %save figure
                saveas(gcf,[path,'plots\headingOffsetVsBumpParameters.png']);
                
            else
                
                figure('Position',[100 100 1000 800]),
                
                subplot(2,2,1)
                binscatter(offset_variability,smooth(allBumpMag))
                correlations = corrcoef(offset_variability,smooth(allBumpMag));
                ylabel('Bump magnitude (max-min)');
                xlabel('Bar offset variability');
                text(1,2.5,['Corr ',num2str(round(correlations(2,1),2))]);
                title('Non-smoothed data');
                
                subplot(2,2,2)
                binscatter(offset_variability,smooth(allHalfWidth))
                correlations = corrcoef(offset_variability,smooth(allHalfWidth));
                ylabel('Bump half width');
                xlabel('Bar offset variability');
                text(1,6,['Corr ',num2str(round(correlations(2,1),2))]);
                title('Non-smoothed data');
                
                %with smoothed offset variability
                subplot(2,2,3)
                binscatter(smoothed_bar_offset_variability,smooth(allBumpMag))
                correlations = corrcoef(smoothed_bar_offset_variability,smooth(allBumpMag));
                ylabel('Bump magnitude (max-min)');
                xlabel('Bar offset variability');
                text(0.8,2.5,['Corr ',num2str(round(correlations(2,1),2))]);
                title('Smoothed data');
                
                subplot(2,2,4)
                binscatter(smoothed_bar_offset_variability,smooth(allHalfWidth))
                correlations = corrcoef(smoothed_bar_offset_variability,smooth(allHalfWidth));
                ylabel('Bump half width');
                xlabel('Bar offset variability');
                text(0.8,6,['Corr ',num2str(round(correlations(2,1),2))]);
                title('Smoothed data');
                
                %save figure
                saveas(gcf,[path,'plots\barOffsetVsBumpParameters.png']);
                
            end
            
            %% Repeat binning the data
            
            if type_of_fly == 1
                
                figure,
                
                nbins = 10;
                max_bin = 1;
                min_bin = 0.2;
                binWidth = max_bin/nbins;
                Bins = [0.2:binWidth:max_bin];
                
                %getting binned means
                for bin = 1:length(Bins)-1
                    meanBin(bin) = median(allBumpMag((heading_offset_variability > Bins(bin)) & (heading_offset_variability < Bins(bin+1))));
                    stdBin(bin) = std(allBumpMag((heading_offset_variability > Bins(bin)) & (heading_offset_variability < Bins(bin+1))));
                    errBin(bin) = stdBin(bin)./sqrt(length(allBumpMag((heading_offset_variability > Bins(bin)) & (heading_offset_variability < Bins(bin+1)))));
                end
                
                %create axes for plot
                mvtAxes = Bins - binWidth;
                mvtAxes = mvtAxes(2:end);
                mvtAxes(end) = mvtAxes(end-1)+binWidth;
                
                %Plot
                subplot(1,2,1)
                boundedline(mvtAxes,meanBin,errBin)
                ylabel({'Bump magnitude';'(amplitude of Fourier component)'}); xlabel('Heading offset variability');
                ylim([0 2.5]);
                title('Bump magnitude');
                
                %getting binned means
                for bin = 1:length(Bins)-1
                    meanBin(bin) = median(allHalfWidth((heading_offset_variability > Bins(bin)) & (heading_offset_variability < Bins(bin+1))));
                    stdBin(bin) = std(allHalfWidth((heading_offset_variability > Bins(bin)) & (heading_offset_variability < Bins(bin+1))));
                    errBin(bin) = stdBin(bin)./sqrt(length(allHalfWidth((heading_offset_variability > Bins(bin)) & (heading_offset_variability < Bins(bin+1)))));
                end
                
                %Plot
                subplot(1,2,2)
                boundedline(mvtAxes,meanBin,errBin)
                ylabel('Bump half width'); xlabel('Heading offset variability');
                ylim([1.5 3.5]);
                title('Bump half width');
                
                %save figure
                saveas(gcf,[path,'plots\headingOffsetVsMedianBumpParameters.png']);
                
            else
                
                figure,
                
                nbins = 10;
                max_bin = 1;
                min_bin = 0.2;
                binWidth = max_bin/nbins;
                Bins = [0.2:binWidth:max_bin];
                
                %getting binned means
                for bin = 1:length(Bins)-1
                    meanBin(bin) = median(allBumpMag((offset_variability > Bins(bin)) & (offset_variability < Bins(bin+1))));
                    stdBin(bin) = std(allBumpMag((offset_variability > Bins(bin)) & (offset_variability < Bins(bin+1))));
                    errBin(bin) = stdBin(bin)./sqrt(length(allBumpMag((offset_variability > Bins(bin)) & (offset_variability < Bins(bin+1)))));
                end
                
                %create axes for plot
                mvtAxes = Bins - binWidth;
                mvtAxes = mvtAxes(2:end);
                mvtAxes(end) = mvtAxes(end-1)+binWidth;
                
                %Plot
                subplot(1,2,1)
                boundedline(mvtAxes,meanBin,errBin)
                ylabel({'Bump magnitude';'(amplitude of Fourier component)'}); xlabel('Bar offset variability');
                ylim([0 2]);
                title('Bump magnitude');
                
                %getting binned means
                for bin = 1:length(Bins)-1
                    meanBin(bin) = median(allHalfWidth((offset_variability > Bins(bin)) & (offset_variability < Bins(bin+1))));
                    stdBin(bin) = std(allHalfWidth((offset_variability > Bins(bin)) & (offset_variability < Bins(bin+1))));
                    errBin(bin) = stdBin(bin)./sqrt(length(allHalfWidth((offset_variability > Bins(bin)) & (offset_variability < Bins(bin+1)))));
                end
                
                %Plot
                subplot(1,2,2)
                boundedline(mvtAxes,meanBin,errBin)
                ylabel('Bump half width'); xlabel('Bar offset variability');
                ylim([1.5 3.5]);
                title('Bump half width');
                
                %save figure
                saveas(gcf,[path,'plots\barOffsetVsMedianBumpParameters.png']);
            end
            
            %% Repeat only for the inverted gain portion
            
            BumpMagIG = allBumpMag(gain_changes(1):gain_changes(2));
            HalfWidthIG = allHalfWidth(gain_changes(1):gain_changes(2));
            heading_offset_variabilityIG = heading_offset_variability(gain_changes(1):gain_changes(2));
            offset_variabilityIG = offset_variability(gain_changes(1):gain_changes(2));
            
            if type_of_fly == 1
                
                figure,
                
                nbins = 10;
                max_bin = 1;
                min_bin = 0.2;
                binWidth = max_bin/nbins;
                Bins = [0.2:binWidth:max_bin];
                
                %getting binned means
                for bin = 1:length(Bins)-1
                    meanBin(bin) = median(BumpMagIG((heading_offset_variabilityIG > Bins(bin)) & (heading_offset_variabilityIG < Bins(bin+1))));
                    stdBin(bin) = std(BumpMagIG((heading_offset_variabilityIG > Bins(bin)) & (heading_offset_variabilityIG < Bins(bin+1))));
                    errBin(bin) = stdBin(bin)./sqrt(length(BumpMagIG((heading_offset_variabilityIG > Bins(bin)) & (heading_offset_variabilityIG < Bins(bin+1)))));
                end
                
                %create axes for plot
                mvtAxes = Bins - binWidth;
                mvtAxes = mvtAxes(2:end);
                mvtAxes(end) = mvtAxes(end-1)+binWidth;
                
                %Plot
                subplot(1,2,1)
                boundedline(mvtAxes,meanBin,errBin)
                ylabel('Bump magnitude (max-min)'); xlabel('Heading offset variability');
                ylim([0 3]);
                title('Bump magnitude');
                
                %getting binned means
                for bin = 1:length(Bins)-1
                    meanBin(bin) = median(HalfWidthIG((heading_offset_variabilityIG > Bins(bin)) & (heading_offset_variabilityIG < Bins(bin+1))));
                    stdBin(bin) = std(HalfWidthIG((heading_offset_variabilityIG > Bins(bin)) & (heading_offset_variabilityIG < Bins(bin+1))));
                    errBin(bin) = stdBin(bin)./sqrt(length(HalfWidthIG((heading_offset_variabilityIG > Bins(bin)) & (heading_offset_variabilityIG < Bins(bin+1)))));
                end
                
                %Plot
                subplot(1,2,2)
                boundedline(mvtAxes,meanBin,errBin)
                ylabel('Bump half width'); xlabel('Heading offset variability');
                ylim([1.5 3.5]);
                title('Bump half width');
                
                %save figure
                saveas(gcf,[path,'plots\headingOffsetVsMedianBumpParametersIG.png']);
                
            else
                
                figure,
                
                nbins = 10;
                max_bin = 1;
                min_bin = 0.2;
                binWidth = max_bin/nbins;
                Bins = [0.2:binWidth:max_bin];
                
                %getting binned means
                for bin = 1:length(Bins)-1
                    meanBin(bin) = median(BumpMagIG((offset_variabilityIG > Bins(bin)) & (offset_variabilityIG < Bins(bin+1))));
                    stdBin(bin) = std(BumpMagIG((offset_variabilityIG > Bins(bin)) & (offset_variabilityIG < Bins(bin+1))));
                    errBin(bin) = stdBin(bin)./sqrt(length(BumpMagIG((offset_variabilityIG > Bins(bin)) & (offset_variabilityIG < Bins(bin+1)))));
                end
                
                %create axes for plot
                mvtAxes = Bins - binWidth;
                mvtAxes = mvtAxes(2:end);
                mvtAxes(end) = mvtAxes(end-1)+binWidth;
                
                %Plot
                subplot(1,2,1)
                boundedline(mvtAxes,meanBin,errBin)
                ylabel({'Bump magnitude';'(amplitude of Fourier component)'}); xlabel('Bar offset variability');
                ylim([0 3]);
                title('Bump magnitude');
                
                %getting binned means
                for bin = 1:length(Bins)-1
                    meanBin(bin) = median(HalfWidthIG((offset_variabilityIG > Bins(bin)) & (offset_variabilityIG < Bins(bin+1))));
                    stdBin(bin) = std(HalfWidthIG((offset_variabilityIG > Bins(bin)) & (offset_variabilityIG < Bins(bin+1))));
                    errBin(bin) = stdBin(bin)./sqrt(length(HalfWidthIG((offset_variabilityIG > Bins(bin)) & (offset_variabilityIG < Bins(bin+1)))));
                end
                
                %Plot
                subplot(1,2,2)
                boundedline(mvtAxes,meanBin,errBin)
                ylabel('Bump half width'); xlabel('Bar offset variability');
                ylim([1.5 3.5]);
                title('Bump half width');
                
                %save figure
                saveas(gcf,[path,'plots\barOffsetVsMedianBumpParametersIG.png']);
            end
            
                        
            %% Heading stability
            
            heading_variability = matlab.tall.movingWindow(fcn,50,deg2rad(heading));
            smoothed_heading_variability = smooth(heading_variability,150,'rloess');
            
            
            figure('Position',[100 100 1600 400]),
            subplot(3,1,1)
            plot(x_out_heading,heading_to_plot,'k')
            hold on
            for change = 1:length(gain_changes)
                line([data.time(gain_changes(change)) data.time(gain_changes(change))], [-180 180], 'LineWidth', 2, 'color', [0,0.5,0]);
            end
            ylim([-180 180]);
            title('Heading');
            
            subplot(3,1,2)
            plot(data.time,smooth(heading_variability),'k')
            hold on
            for change = 1:length(gain_changes)
                line([data.time(gain_changes(change)) data.time(gain_changes(change))], [0 2], 'LineWidth', 2, 'color', [0,0.5,0]);
            end
            title('Heading variability');
            
            subplot(3,1,3)
            plot(data.time,smoothed_heading_variability,'k')
            hold on
            for change = 1:length(gain_changes)
                line([data.time(gain_changes(change)) data.time(gain_changes(change))], [0 2], 'LineWidth', 2, 'color', [0,0.5,0]);
            end
            title('Smoothed heading variability');
            
            %save figure
            saveas(gcf,[path,'plots\HeadingVariability.png']);
            
            
            %% plot relationship between the offset and the heading variability
            
            figure('Position',[100 100 1000 800]),
            
            subplot(2,2,1)
            binscatter(heading_variability,offset_variability)
            correlations = corrcoef(heading_variability,offset_variability);
            ylabel('Bar offset variability');
            xlabel('Heading variability');
            text(1,1.2,['Corr ',num2str(round(correlations(2,1),2))]);
            title('Non-smoothed data');
            
            subplot(2,2,2)
            binscatter(heading_variability,heading_offset_variability)
            correlations = corrcoef(heading_variability,heading_offset_variability);
            ylabel('Heading offset variability');
            xlabel('Heading variability');
            text(1,1.2,['Corr ',num2str(round(correlations(2,1),2))]);
            title('Non-smoothed data');
            
            %with smoothed offset variability
            subplot(2,2,3)
            binscatter(smoothed_heading_variability,smoothed_bar_offset_variability)
            correlations = corrcoef(smoothed_heading_variability,smoothed_bar_offset_variability);
            ylabel('Bar offset variability');
            xlabel('Heading variability');
            text(1,1.2,['Corr ',num2str(round(correlations(2,1),2))]);
            title('Smoothed data');
            
            subplot(2,2,4)
            binscatter(smoothed_heading_variability,smoothed_heading_offset_variability)
            correlations = corrcoef(smoothed_heading_variability,smoothed_heading_offset_variability);
            ylabel('Heading offset variability');
            xlabel('Heading variability');
            text(1,1.2,['Corr ',num2str(round(correlations(2,1),2))]);
            title('Smoothed data');
            
            %save figure
            saveas(gcf,[path,'plots\headingVarVsOffsetVar.png']);
            
            
            %% Models of bump parameters
            %we will model the two bump parameters as a function of total movement and
            %offset variability
            
            %zscore the total movement data
            %fill missing data in the total movement vector
            data.total_mvt_ds = fillmissing(data.total_mvt_ds,'linear');
            zscored_total_mvt = zscore(data.total_mvt_ds);
            
            %crate table with the model's variables
            modelTable = table(offset_variability,heading_offset_variability,smoothed_bar_offset_variability,smoothed_heading_offset_variability,data.total_mvt_ds(1:length(allBumpMag))',zscored_total_mvt(1:length(allBumpMag))',data.time(1:length(allBumpMag)),allBumpMag',allHalfWidth',allZBumpMag',allZHalfWidth',heading_variability,'VariableNames',{'BarOffsetVariability','HeadingOffsetVariability','SmoothBarOffsetVariability','SmoothHeadingOffsetVariability','TotalMovement','ZscoredTotalMovement','Time','BumpMagnitude','BumpWidth','ZBumpMagnitude','ZBumpWidth','HeadingVariability'});
            
            if type_of_fly == 1
                mdl_BM = fitlm(modelTable,'BumpMagnitude~HeadingOffsetVariability+ZscoredTotalMovement+Time')
            else
                mdl_BM = fitlm(modelTable,'BumpMagnitude~BarOffsetVariability+ZscoredTotalMovement+Time')
            end
            
            
            %% for bump width at half max
            
            if type_of_fly == 1
                mdl_HW = fitlm(modelTable,'BumpWidth~HeadingOffsetVariability+ZscoredTotalMovement+Time')
            else
                mdl_HW = fitlm(modelTable,'BumpWidth~BarOffsetVariability+ZscoredTotalMovement+Time')
            end
            

            
            %% Save relevant data
            
            save([path,'\gain_change_data.mat'],'modelTable','type_of_fly');
            
            close all;
            
        end
    end
end
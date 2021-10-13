%code to run the offset verifier block on all flies of the dataset

function run_offset_verifier_block_analysis()
%% Load data

%Get the path for each fly
parentDir = 'Z:\Wilson Lab\Mel\Experiments\Uncertainty\Exp25\data\Experimental\two_ND_filters_3_contrasts';
folderNames = dir(parentDir);

for content = 1:length(folderNames)
    if contains(folderNames(content).name,'60D05')
        flyData{content} = [folderNames(content).folder,'\',folderNames(content).name];
    end
end

%remove empty cells
data_dirs = flyData(~cellfun(@isempty,flyData));

%for each fly, load the data from all the open loop sessions, analyze and
%plot
for fly = 1:length(data_dirs)
    
    load([data_dirs{fly},'\sessions_info.mat']);
    sid = sessions_info.offset_verifier;
    
    %get the contents of the fly folder
    fly_files = dir([data_dirs{fly},'\analysis']);
    %determine which content belongs to the sid we extracted
    
    %determine which content belongs to the sid we extracted
    for file = 1:length(fly_files)
        if (contains(fly_files(file).name,['sid_',num2str(sid),'_']) & contains(fly_files(file).name,'continuous'))
            %load the data
            load([fly_files(file).folder,'\',fly_files(file).name])
                      
            %% Make directory to save plots
            
            path = [fly_files(file).folder,'\'];
            %Move to the analysis folder
            cd(path)
            %List the contents
            contents = dir();
            %if there isn't a 'plots' folder already, create one
            if (contains([contents.name],'continuous_plots') == 0)
                mkdir(path,'continuous_plots');
            end
            %List the contents of the 'plots' folder
            cd([path,'continuous_plots\'])
            
            %% Plot heatmap with bar position, phase and offset
            
            figure('Position',[200 200 1000 600]),
            subplot(3,6,[1 5])
            %Plot heatmap of EPG activity
            dff_matrix = continuous_data.dff_matrix';
            imagesc(flip(dff_matrix))
            colormap(flipud(gray))
            set(gca,'ytick',[]);
            set(gca,'yticklabel',{[]});
            title('EPG activity in the PB');
            set(gca,'xtick',[]);
            set(gca,'xticklabel',{[]});
            
            subplot(3,6,[7 11])
            %Get bar position to plot
            bar_position = wrapTo180(continuous_data.panel_angle);
            [x_out_bar,bar_position_to_plot] = removeWrappedLines(continuous_data.time,bar_position);
            plot(x_out_bar,bar_position_to_plot,'LineWidth',1.5,'color',[0.2 0.6 0.7])
            hold on
            %Get EPG phase to plot
            phase = wrapTo180(rad2deg(continuous_data.bump_pos'));
            [x_out_phase,phase_to_plot] = removeWrappedLines(continuous_data.time,phase);
            plot(x_out_phase,phase_to_plot,'LineWidth',1.5,'color',[0.9 0.3 0.4])
            xlim([0 x_out_phase(end-1)]);
            ylim([-180 180]);
            axP = get(gca,'Position');
            legend('Bar position','EPG phase','Location','EastOutside')
            set(gca, 'Position', axP);
            ylabel('Deg');
            title('Bar and bump position');
            set(gca,'xticklabel',{[]})
            
            subplot(3,6,[13 17])
            %Get offset to plot
            offset = wrapTo180(continuous_data.offset);
            [x_out_offset,offset_to_plot] = removeWrappedLines(continuous_data.time,offset);
            plot(x_out_offset,offset_to_plot,'k','LineWidth',1.5)
            xlim([0 x_out_offset(end-1)]);
            ylim([-180 180]);
            xlabel('Time (sec)');
            ylabel('Deg');
            title('Offset');
            
            subplot(3,6,18)
            polarhistogram(deg2rad(offset),'FaceColor','k')
            set(gca,'ThetaZeroLocation','top',...
                'ThetaDir','counterclockwise');
            title('Offset distribution');
            rticklabels({}); %remove rticklabels
            
            %Save figure
            saveas(gcf,[path,'continuous_plots\Offset_verifyier_block_heatmap.png']);
            
            %% Get reference offset
            
            visual_offset = circ_dist(continuous_data.bump_pos',deg2rad(continuous_data.panel_angle));
            mean_verifying_offset = rad2deg(circ_mean(visual_offset,[],1));
            
            %save
            save([path,'continuous_verifying_offset.mat'],'mean_verifying_offset');
            
            
            %% Analyze relationship between yaw speed data and BM
            
            nbins = 20;
            yaw_speed = abs(continuous_data.vel_yaw_ds);
            maxBinYS = min([nanmean(yaw_speed)+3*nanstd(yaw_speed),max(yaw_speed)]);
            binWidthYS = maxBinYS/nbins;
            YSBins = [0:binWidthYS:maxBinYS];
            
            %getting binned means
            for bin = 1:length(YSBins)-1
                meanBin(bin) = mean(continuous_data.bump_magnitude((yaw_speed(1:length(continuous_data.bump_magnitude)) > YSBins(bin)) & (yaw_speed(1:length(continuous_data.bump_magnitude)) < YSBins(bin+1))));
            end
            
            %create axes for plot
            YSAxes = YSBins - binWidthYS;
            YSAxes = YSAxes(2:end);
            YSAxes(end) = YSAxes(end-1)+binWidthYS;
            
            %Plot
            figure('Position',[200 200 1400 600]),
            %Plot bump magnitude in time
            subplot(2,4,[1 3])
            plot(continuous_data.time(1:length(continuous_data.bump_magnitude)),continuous_data.bump_magnitude, 'k')
            ylabel('Bump magnitude');
            ylim([0 5]);
            set(gca,'xticklabel',{[]});
            
            %Plot total movement
            subplot(2,4,[5 7])
            plot(continuous_data.time(1:length(continuous_data.bump_magnitude)),yaw_speed(1:length(continuous_data.bump_magnitude)),'k')
            xlabel('Time (sec)');
            ylabel('Yaw speed (deg/s)');
            
            %Plot relationship between both parameters
            subplot(2,4,[4,8]);
            plot(YSAxes,meanBin,'-ko')
            ylabel('Mean bump magnitude'); xlabel('Yaw speed (deg/s)');
            ylim([0 (max(meanBin)+0.5)]);
            
            
            %% Bin for vel data
            
            for_vel = abs(continuous_data.vel_for_ds');
            maxBinFV = min([nanmean(for_vel)+3*nanstd(for_vel),max(for_vel)]);
            binWidthFV = maxBinFV/nbins;
            FVBins = [0:binWidthFV:maxBinFV];
            
            %getting binned means
            for bin = 1:length(FVBins)-1
                meanBin(bin) = mean(continuous_data.bump_magnitude((for_vel(1:length(continuous_data.bump_magnitude)) > FVBins(bin)) & (for_vel(1:length(continuous_data.bump_magnitude)) < FVBins(bin+1))));
            end
            
            %create axes for plot
            fvAxes = FVBins - binWidthFV;
            fvAxes = fvAxes(2:end);
            fvAxes(end) = fvAxes(end-1)+binWidthFV;
            
            %Plot
            figure('Position',[200 200 1400 600]),
            %Plot bump magnitude in time
            subplot(2,4,[1 3])
            plot(continuous_data.time(1:length(continuous_data.bump_magnitude)),continuous_data.bump_magnitude, 'k')
            ylabel('Bump magnitude');
            ylim([0 5]);
            set(gca,'xticklabel',{[]});
            
            %Plot total movement
            subplot(2,4,[5 7])
            plot(continuous_data.time(1:length(continuous_data.bump_magnitude)),for_vel(1:length(continuous_data.bump_magnitude)),'k')
            xlabel('Time (sec)');
            ylabel('Forward velocity (mm/s)');
            
            %Plot relationship between both parameters
            subplot(2,4,[4,8]);
            plot(fvAxes,meanBin,'-ko')
            ylabel('Mean bump magnitude'); xlabel('Forward velocity (mm/s)');
            ylim([0 (max(meanBin)+0.5)]);
            
            %% Bin both parameters to build heatmap
            
            %getting binned means
            for ys_bin = 1:length(YSBins)-1
                for fv_bin = 1:length(FVBins)-1
                    doubleBin(ys_bin,fv_bin) = nanmean(continuous_data.bump_magnitude((yaw_speed(1:length(continuous_data.bump_magnitude)) > YSBins(ys_bin)) & (yaw_speed(1:length(continuous_data.bump_magnitude)) < YSBins(ys_bin+1)) & (for_vel(1:length(continuous_data.bump_magnitude)) > FVBins(fv_bin)) & (for_vel(1:length(continuous_data.bump_magnitude)) < FVBins(fv_bin+1))));
                end
            end
            
            %flip the data such that high forward velocity values are at the top
            binned_data = flip(doubleBin);
            
            figure,
            imagesc(binned_data)
            xticks([1:4:20])
            set(gca, 'XTickLabel', round(FVBins(1:4:20)))
            xlabel('Forward velocity (mm/s)','fontsize',12,'fontweight','bold');
            ylabel('Yaw speed (deg/sec)','fontsize',12,'fontweight','bold');
            yticks([1:4:20])
            set(gca, 'YTickLabel', round(YSBins(20:-4:1)))
            c = colorbar;
            ylabel(c, 'Bump magnitude')
            
            
            %Save figure
            saveas(gcf,[path,'continuous_plots\bm_heatmap_offset_verifying_block.png']);
            
            bump_magnitude = continuous_data.bump_magnitude;
            save([path,'\continuous_vel_bm_data_ov.mat'],'for_vel','yaw_speed','bump_magnitude');
            
            
        end
    end
    
    close all;
end
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
        if contains(fly_files(file).name,['sid_',num2str(sid),'_'])
            %load the data
            load([fly_files(file).folder,'\',fly_files(file).name])
            
            
            %% Make directory to save plots
            
            path = [fly_files(file).folder,'\'];
            %Move to the analysis folder
            cd(path)
            %List the contents
            contents = dir();
            %if there isn't a 'plots' folder already, create one
            if (contains([contents.name],'plots') == 0)
                mkdir(path,'plots');
            end
            %List the contents of the 'plots' folder
            cd([path,'plots\'])
            %% Plot heatmap with bar position, phase and offset
            
            figure('Position',[200 200 1000 600]),
            subplot(3,6,[1 5])
            %Plot heatmap of EPG activity
            imagesc(data.dff_matrix)
            colormap(flipud(gray))
            ylabel('PB glomerulus');
            yticks(1:2:16);
            yticklabels({'8L','6L','4L','2L','1R','3R','5R','7R'});
            title('EPG activity in the PB');
            set(gca,'xtick',[]);
            set(gca,'xticklabel',{[]});
            
            subplot(3,6,[7 11])
            %Get bar position to plot
            bar_position = wrapTo180(data.panel_angle);
            [x_out_bar,bar_position_to_plot] = removeWrappedLines(data.time,bar_position);
            plot(x_out_bar,bar_position_to_plot,'LineWidth',1.5,'color',[0.2 0.6 0.7])
            hold on
            %Get EPG phase to plot
            phase = wrapTo180(rad2deg(data.dff_pva));
            [x_out_phase,phase_to_plot] = removeWrappedLines(data.time,phase);
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
            offset = wrapTo180(data.offset);
            [x_out_offset,offset_to_plot] = removeWrappedLines(data.time,offset);
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
            saveas(gcf,[path,'plots\Offset_verifyier_block_heatmap.png']);
            
            %% Get reference offset
            
            visual_offset = circ_dist(data.dff_pva,deg2rad(data.panel_angle));
            mean_verifying_offset = rad2deg(circ_mean(visual_offset,[],1));
            
            %save
            save([path,'\verifying_offset.mat'],'mean_verifying_offset');
            
                        
            %% Analyze the relationship between fly velocity and BM
  
            nbins = 10;
            yaw_speed = abs(data.vel_yaw_ds);
            maxBinYS = min([max(yaw_speed), nanmean(yaw_speed)+3*nanstd(yaw_speed)]);
            binWidthYS = maxBinYS/nbins;
            YSBins10 = [0:binWidthYS:maxBinYS];
            for_vel = abs(data.vel_for_ds');
            maxBinFV = min([max(for_vel), nanmean(for_vel)+3*nanstd(for_vel)]);
            binWidthFV = maxBinFV/nbins;
            FVBins10 = [0:binWidthFV:maxBinFV];
            
            %getting binned means
            for ys_bin = 1:length(YSBins10)-1
                for fv_bin = 1:length(FVBins10)-1
                    doubleBin10(ys_bin,fv_bin) = nanmean(data.bump_magnitude((yaw_speed(1:length(data.bump_magnitude)) > YSBins10(ys_bin)) & (yaw_speed(1:length(data.bump_magnitude)) < YSBins10(ys_bin+1)) & (for_vel(1:length(data.bump_magnitude)) > FVBins10(fv_bin)) & (for_vel(1:length(data.bump_magnitude)) < FVBins10(fv_bin+1))));
                end
            end
            
            %flip the data such that high forward velocity values are at the top
            binned_data10 = flip(doubleBin10);
            
            figure,
            imagesc(binned_data10)
            xticks([1:2:10])
            set(gca, 'XTickLabel', round(FVBins10(1:2:10)))
            xlabel('Forward velocity (mm/s)','fontsize',12,'fontweight','bold');
            ylabel('Yaw speed (deg/sec)','fontsize',12,'fontweight','bold');
            yticks([1:2:10])
            set(gca, 'YTickLabel', round(YSBins10(10:-2:1)))
            c = colorbar;
            ylabel(c, 'Bump magnitude')
            
            %Save figure
            saveas(gcf,[fly_files(file).folder,'\plots\bm_heatmap_10bins_ov.png']);   
            
            %% Repeat with more bins
            
            nbins = 20;
            binWidthYS = maxBinYS/nbins;
            YSBins = [0:binWidthYS:maxBinYS];
            binWidthFV = maxBinFV/nbins;
            FVBins = [0:binWidthFV:maxBinFV];
            
            %getting binned means
            for ys_bin = 1:length(YSBins)-1
                for fv_bin = 1:length(FVBins)-1
                    doubleBin(ys_bin,fv_bin) = nanmean(data.bump_magnitude((yaw_speed(1:length(data.bump_magnitude)) > YSBins(ys_bin)) & (yaw_speed(1:length(data.bump_magnitude)) < YSBins(ys_bin+1)) & (for_vel(1:length(data.bump_magnitude)) > FVBins(fv_bin)) & (for_vel(1:length(data.bump_magnitude)) < FVBins(fv_bin+1))));
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
            saveas(gcf,[fly_files(file).folder,'\plots\bm_heatmap_ov.png']);   
            
            %Save velocity and bump magnitude data
            bump_magnitude = data.bump_magnitude;
            save([fly_files(file).folder,'\vel_bm_data_ov.mat'],'for_vel','yaw_speed','bump_magnitude');
            
            close all
            
        end
    end
    
    close all;
end
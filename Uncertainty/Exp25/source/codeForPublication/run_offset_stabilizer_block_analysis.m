
function run_offset_stabilizer_block_analysis()

%This function plots the data for the offset stabilizer block for all the flies in the
%dataset (heatmap of EPG activity, phase and bar position in time, offset
%in time, and offset distribution)



%Get the path to the data for each fly
parentDir = 'Z:\Wilson Lab\Mel\Experiments\Uncertainty\Exp25\data\Experimental\two_ND_filters_3_contrasts';
folderNames = dir(parentDir);

for content = 1:length(folderNames)
    if contains(folderNames(content).name,'60D05')
        flyData{content} = [folderNames(content).folder,'\',folderNames(content).name];
    end
end
%Remove empty cells
data_dirs = flyData(~cellfun(@isempty,flyData));


%Determine the session id for this trial type from the session_info.mat
%file, load the data and plot
for fly = 1:length(data_dirs)
    
    load([data_dirs{fly},'\sessions_info.mat']);
    sid = sessions_info.offset_stabilizer;
    
    %Get the contents of the fly folder
    fly_files = dir([data_dirs{fly},'\analysis']);
    
    %Determine which content belongs to the sid we extracted
    for file = 1:length(fly_files)
        if contains(fly_files(file).name,['sid_',num2str(sid),'_'])
            %Load the data
            load([fly_files(file).folder,'\',fly_files(file).name])
            
            %% Plot heatmap with bar position, phase and offset
            
            figure('Position',[200 200 1000 600]),
            subplot(3,6,[1 5])
            %Plot heatmap of EPG activity
            imagesc(data.dff_matrix)
            colormap(flipud(gray))
            ylabel('PB glomerulus','fontweight','bold','fontsize',10);
            yticks(1:2:16);
            yticklabels({'8L','6L','4L','2L','1R','3R','5R','7R'});
            title('EPG activity in the PB','fontweight','bold','fontsize',12);
            set(gca,'xtick',[]);
            set(gca,'xticklabel',{[]});
            
            subplot(3,6,[7 11])
            %Get bar position to plot
            bar_position = wrapTo180(data.panel_angle);
            %Remove lines that wrap around using auxiliary function
            [x_bar_pos,bar_position_to_plot] = removeWrappedLines(data.time,bar_position);
            plot(x_bar_pos,bar_position_to_plot,'LineWidth',1.5,'color',[0.2 0.6 0.7])
            hold on
            %Get EPG phase to plot
            phase = wrapTo180(rad2deg(data.dff_pva));
            %Remove lines that wrap around using auxiliary function
            [x_out_phase,phase_to_plot] = removeWrappedLines(data.time,phase);
            plot(x_out_phase,phase_to_plot,'LineWidth',1.5,'color',[0.9 0.3 0.4])
            if ~isnan(x_out_phase(end))
                xlim([0 x_out_phase(end)]);
            else
                %there's one fly in which the last frame is a NaN, so use the previous one in that case
                xlim([0 x_out_phase(end-1)]); 
            end
            ylim([-180 180]);
            legend('Bar position','EPG phase','Location','EastOutside')
            %Set legend outside of the plot so that it doesn't occlude the traces
            axP = get(gca,'Position');
            set(gca, 'Position', axP);
            ylabel('Deg','fontweight','bold','fontsize',10);
            title('Bar and bump position','fontweight','bold','fontsize',12);
            set(gca,'xticklabel',{[]})
            
            subplot(3,6,[13 17])
            %Get offset to plot
            offset = wrapTo180(data.offset);
            %Remove lines that wrap around using auxiliary function
            [x_out_offset,offset_to_plot] = removeWrappedLines(data.time,offset);
            plot(x_out_offset,offset_to_plot,'k','LineWidth',1.5)
            xlim([0 x_out_offset(end)]);
            ylim([-180 180]);
            xlabel('Time (sec)','fontweight','bold','fontsize',10);
            ylabel('Deg','fontweight','bold','fontsize',10);
            title('Offset','fontweight','bold','fontsize',12);
            
            %Plot circular offset distribution
            subplot(3,6,18)
            polarhistogram(deg2rad(offset),'FaceColor','k')
            set(gca,'ThetaZeroLocation','top',...
                'ThetaDir','counterclockwise');
            title('Offset distribution','fontweight','bold','fontsize',12);
            rticklabels({}); %remove rticklabels
            
            %Save figure
            saveas(gcf,[fly_files(file).folder,'\plots\Offset_stabilizer_block_heatmap.png']);
            
            
            %% Analyze the relationship between fly velocity and BM
  
            nbins = 10;
            yaw_speed = abs(data.vel_yaw_ds);
            maxBinYS = max(yaw_speed);
            binWidthYS = maxBinYS/nbins;
            YSBins10 = [0:binWidthYS:maxBinYS];
            for_vel = abs(data.vel_for_ds');
            maxBinFV = max(for_vel);
            binWidthFV = maxBinFV/nbins;
            FVBins10 = [0:binWidthFV:maxBinFV];
            
            %getting binned means
            for ys_bin = 1:length(YSBins10)-1
                for fv_bin = 1:length(FVBins10)-1
                    doubleBin10(ys_bin,fv_bin) = mean(data.bump_magnitude((yaw_speed(1:length(data.bump_magnitude)) > YSBins10(ys_bin)) & (yaw_speed(1:length(data.bump_magnitude)) < YSBins10(ys_bin+1)) & (for_vel(1:length(data.bump_magnitude)) > FVBins10(fv_bin)) & (for_vel(1:length(data.bump_magnitude)) < FVBins10(fv_bin+1))));
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
            saveas(gcf,[fly_files(file).folder,'\plots\bm_heatmap_10bins.png']);   
            
            %% Repeat with more bins
            
            nbins = 20;
            binWidthYS = maxBinYS/nbins;
            YSBins = [0:binWidthYS:maxBinYS];
            binWidthFV = maxBinFV/nbins;
            FVBins = [0:binWidthFV:maxBinFV];
            
            %getting binned means
            for ys_bin = 1:length(YSBins)-1
                for fv_bin = 1:length(FVBins)-1
                    doubleBin(ys_bin,fv_bin) = mean(data.bump_magnitude((yaw_speed(1:length(data.bump_magnitude)) > YSBins(ys_bin)) & (yaw_speed(1:length(data.bump_magnitude)) < YSBins(ys_bin+1)) & (for_vel(1:length(data.bump_magnitude)) > FVBins(fv_bin)) & (for_vel(1:length(data.bump_magnitude)) < FVBins(fv_bin+1))));
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
            saveas(gcf,[fly_files(file).folder,'\plots\bm_heatmap.png']);   
            
            %Save velocity and bump magnitude data
            bump_magnitude = data.bump_magnitude;
            save([fly_files(file).folder,'\vel_bm_data.mat'],'for_vel','yaw_speed','bump_magnitude');
            
            close all
            
        end
                
    end
end
end

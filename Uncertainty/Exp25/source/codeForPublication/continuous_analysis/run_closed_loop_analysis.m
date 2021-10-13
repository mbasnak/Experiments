
function run_closed_loop_analysis()

%run the closed loop analysis for all the flies in the
%dataset

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

%determine the session id for this trial type from the session_info.mat
%file, load the data and plot
for fly = 1:length(data_dirs)
    
    load([data_dirs{fly},'\sessions_info.mat']);
    sid = sessions_info.closed_loop;
    
    %get the contents of the fly folder
    fly_files = dir([data_dirs{fly},'\analysis']);
    %determine which content belongs to the sid we extracted
    
    %determine which content belongs to the sid we extracted
    for file = 1:length(fly_files)
        if (contains(fly_files(file).name,['sid_',num2str(sid),'_']) & contains(fly_files(file).name,'continuous'))
            %load the data
            fileName = fly_files(file).name;
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
            
            
            
            %% Identify changes in stimulus
            
            %Plot the signal from the y dimension of the panels (since we are
            %controlling the change in stimulus through a y panel function)
            figure,
            subplot 121,
            plot(continuous_data.fr_y_ds)
            title('Panels y signal');
            xlabel('Time (downsampled frames)');
            ylabel('Voltage');
            xlim([0 length(continuous_data.fr_y_ds)]);
            
            subplot 122,
            change_y_panels = abs(diff(continuous_data.fr_y_ds));
            plot(change_y_panels)
            %Obtain the frames where the stim changes using the derivative of the
            %signal from the y dimension of the panels
            changeContrast = find(abs(diff(continuous_data.fr_y_ds))>1);
            hold on
            %Add them to the plot
            plot(changeContrast,change_y_panels(changeContrast),'ro')
            title('Change in y panels');
            xlim([0 length(continuous_data.fr_y_ds)]);
            xlabel('Time (downsampled frames)');
            
            
            %% Identify contrast order
            
            %Recover the position function we used to change the pattern
            pos_function = continuous_data.run_obj.function_number;
            
            %Define the order of intensities of the visual stimuli presented according to the function used
            %1 = darkness; 2 = low contrast, 3 = high contrast
            if pos_function == 195
                Intensities = [1,2,1,3,2,3];
            elseif pos_function == 196
                Intensities = [2,1,3,1,2,3];
            else
                Intensities = [3,1,2,1,2,3];
            end
            
            %For the following fly, cut the last intensity, since the trial
            %malfunctioned in the last part (we have ommitted that data part)
            if (contains(path, '20201020_60D05_7f') & ~contains(path,'fly2'))
                Intensities = Intensities(1:5);
            end
            
            %% Plot the activity heatmap, phase and heading positions, and offset in time
            
            %Define the number of columns in the subplots depending on the fly (the one
            %where fictrac malfunctioned will have one less column)
            num_subplots = length(Intensities);
            
            %Plot
            % Plot the heatmap of EPG activity
            figure('Position',[0 0 1800 800]),
            subplot(4,num_subplots,[1 num_subplots])
            dff_matrix = continuous_data.dff_matrix';
            imagesc(flip(dff_matrix))
            colormap(flipud(gray))
            hold on
            %add the changes in stim
            for change = 1:length(changeContrast)
                line([changeContrast(change) changeContrast(change)], [0 size(dff_matrix,1)], 'LineWidth', 3, 'color', [0, 0.5, 0]);
            end
            set(gca,'YTickLabel',[]);
            title('EPG activity in the PB','fontweight','bold','fontsize',12);
            set(gca,'XTickLabel',[]);
            legend('Change in stimulus');
            
            % Plot the heading and the EPG phase
            subplot(4,num_subplots,[num_subplots+1 num_subplots*2])
            %Get heading to plot
            heading = wrapTo180(continuous_data.heading_deg);
            %Remove wrapped lines to plot
            [x_out_heading,heading_to_plot] = removeWrappedLines(continuous_data.time,heading);
            plot(x_out_heading,heading_to_plot,'LineWidth',1.5,'color',[0.2 0.6 0.7])
            hold on
            %Get EPG phase to plot
            %I'm now going to negate the phase, since I'm plotting heading instead of
            %bar position, so the bump moves in the other direction
            phase = wrapTo180(rad2deg(-continuous_data.bump_pos'));
            [x_out_phase,phase_to_plot] = removeWrappedLines(continuous_data.time,phase);
            plot(x_out_phase,phase_to_plot,'color',[0.9 0.3 0.4],'LineWidth',1.5)
            %add the changes in stim
            for change = 1:length(changeContrast)
                line([continuous_data.time(changeContrast(change)) continuous_data.time(changeContrast(change))], [-180 180], 'LineWidth', 3, 'color', [0, 0.5, 0]);
            end
            legend('Fly heading', 'EPG phase');
            title('Bar and bump position','fontweight','bold','fontsize',12);
            ylim([-180, 180]);
            if ~isnan(x_out_phase(end))
                xlim([1,x_out_phase(end)]);
            else
                xlim([1,x_out_phase(end-1)]);
            end
            ylabel('Deg','fontweight','bold','fontsize',10);
            set(gca,'XTickLabel',[]);
            set(gca,'XTick',[]);
            
            % Plot the offset
            subplot(4,num_subplots,[num_subplots*2+1 num_subplots*3])
            %Get offset to plot
            offset = wrapTo180(continuous_data.offset);
            [x_out_offset,offset_to_plot] = removeWrappedLines(continuous_data.time,offset);
            plot(x_out_offset,offset_to_plot,'LineWidth',1.5,'color','k')
            %Add the changes in stim
            for change = 1:length(changeContrast)
                line([continuous_data.time(changeContrast(change)) continuous_data.time(changeContrast(change))], [-180 180], 'LineWidth', 3, 'color', [0, 0.5, 0]);
            end
            ylim([-180 180]);
            ylabel('Deg','fontweight','bold','fontsize',10); xlabel('Time (sec)','fontweight','bold','fontsize',10);
            title('Offset','fontweight','bold','fontsize',12);
            if ~isnan(x_out_offset(end))
                xlim([1 x_out_offset(end)]);
            else
                xlim([1 x_out_offset(end-1)]);
            end
            
            % Polar histograms of offset
            % Color histograms acording to the intensity level of the bar, using the
            % vector Intensities
            color_gradient = {[0,0,0],[0 0 0.6],[ 0.5 0.8 0.9]};
            subplot(4,num_subplots,num_subplots*3+1)
            polarhistogram(deg2rad(offset(1:changeContrast(1))),15,'FaceColor',color_gradient{Intensities(1)})
            set(gca,'ThetaZeroLocation','top',...
                'ThetaDir','counterclockwise');
            if Intensities(1) == 1
                title({'Offset distribution','darkness'});
            elseif Intensities(1) == 2
                title({'Offset distribution','low contrast'});
            else
                title({'Offset distribution','high contrast'});
            end
            
            for contrast = 1:length(changeContrast)-1
                subplot(4,num_subplots,num_subplots*3+1+contrast)
                polarhistogram(deg2rad(offset(changeContrast(contrast):changeContrast(contrast+1))),15,'FaceColor',color_gradient{Intensities(contrast+1)})
                set(gca,'ThetaZeroLocation','top',...
                    'ThetaDir','counterclockwise');
                if Intensities(contrast+1) == 1
                    title({'Offset distribution','darkness'});
                elseif Intensities(contrast+1) == 2
                    title({'Offset distribution','low contrast'});
                else
                    title({'Offset distribution','high contrast'});
                end
            end
            
            %add last contrast for pos function 197, which appears to have one less
            %change in stim
            if pos_function == 197
                subplot(4,num_subplots,num_subplots*4)
                polarhistogram(deg2rad(offset(changeContrast(contrast+1):end)),15,'FaceColor',color_gradient{Intensities(num_subplots)})
                set(gca,'ThetaZeroLocation','top',...
                    'ThetaDir','counterclockwise');
                if Intensities(num_subplots) == 3
                    title({'Offset distribution','high contrast'});
                elseif Intensities(num_subplots) == 2
                    title({'Offset distribution','low contrast'});
                end
            end
            
            %save figure
            saveas(gcf,[path,'continuous_plots\closedLoopHeatmapAndOffset.png']);
            
            %% Set block limits: initial and last frame of each contrast block
            
            blockLimits{1} = [1,changeContrast(1)-1];
            if length(changeContrast) == 5
                for block = 2:5
                    blockLimits{block} = [changeContrast(block-1),changeContrast(block)-1];
                end
                blockLimits{6} = [changeContrast(5),length(heading)];
            elseif length(changeContrast) == 6
                for block = 2:6
                    blockLimits{block} = [changeContrast(block-1),changeContrast(block)-1];
                end
            else
                for block = 2:4
                    blockLimits{block} = [changeContrast(block-1),changeContrast(block)-1];
                end
                blockLimits{5} = [changeContrast(4),length(heading)];
            end
            
            %% Plot offset variation (i.e., compass error)
            
            figure('Position',[200 200 1000 800]),
            subplot(2,1,1)
            for block = 1:length(blockLimits)
                %Compute offset variation
                [~, offset_var(block)] = circ_std(deg2rad(offset(blockLimits{block}(1):blockLimits{block}(2))),[],[],1);
                plot(block,offset_var(block),'ko','MarkerFaceColor',color_gradient{Intensities(block)},'MarkerSize',8)
                hold on
            end
            xlim([0 num_subplots+1]);
            ylim([min(offset_var)-0.5 max(offset_var)+0.5]);
            xticks(1:num_subplots);
            title('Offset variation per block');
            ylabel({'Circular standard deviation','of the offset'});
            xlabel('Block number');
            %Add custom legend with the appropriate colors
            h = zeros(3, 1);
            h(1) = plot(NaN,NaN,'ko','MarkerFaceColor',color_gradient{1});
            h(2) = plot(NaN,NaN,'ko','MarkerFaceColor',color_gradient{2});
            h(3) = plot(NaN,NaN,'ko','MarkerFaceColor',color_gradient{3});
            legend(h, 'Darkness','Low contrast','High Contrast','location','best');
            
            subplot(2,1,2)
            %Create table with contrast level and offset variation values
            for block = 1:length(blockLimits)
                if Intensities(block) == 1
                    contrasts{block} = 'Darkness';
                elseif Intensities(block) == 2
                    contrasts{block} = 'Low contrast';
                else
                    contrasts{block} = 'High contrast';
                end
            end
            summary_data = table(contrasts',offset_var','VariableNames',{'contrast','offset_var'});
            %Get mean offset var by contrast
            mean_data = varfun(@mean,summary_data,'InputVariables','offset_var',...
                'GroupingVariables',{'contrast'});
            %Plot
            %The mean_data sorts contrasts alphabetically, so high contrast appears
            %second, and low contrast third. We will take this into account by creating
            %a 'newOrder' variable
            newOrder = [1,3,2];
            for contrastLevel = 1:3
                plot(contrastLevel,mean_data.mean_offset_var(newOrder(contrastLevel)),'ko','MarkerFaceColor',color_gradient{contrastLevel},'MarkerSize',8)
                hold on
            end
            xlim([0 4]);
            ylim([min(mean_data.mean_offset_var)-0.5 max(mean_data.mean_offset_var)+0.5]);
            xticks(1:3);
            xticklabels({'Darkness', 'Low contrast', 'High contrast'});
            title('Mean offset variation per contrast');
            ylabel({'Circular standard deviation','of the offset'});
            
            %save figure
            saveas(gcf,[path,'continuous_plots\closedLoopOffsetVariation.png']);
            
            %% Get offset value from last bout of high contrast: this will be our 'reference offset' to use during the open loop analysis
            
            %Find limits of last high contrast bout
            high_contrast_bouts = find(Intensities(1:length(blockLimits))==3);
            last_high_contrast = high_contrast_bouts(end);
            
            %This will be the 'visual' offset, so we need to compute it as the circular
            %distance between the stimulus and the phase
            visual_offset = circ_dist(continuous_data.bump_pos',deg2rad(continuous_data.panel_angle));
            visual_offset2 = circ_dist(-continuous_data.bump_pos',deg2rad(-continuous_data.panel_angle));            
            
            %Compute the circular mean of this visual offset in the last high contrast
            %bout
            mean_reference_offset = rad2deg(circ_mean(visual_offset(blockLimits{last_high_contrast}(1):blockLimits{last_high_contrast}(2)),[],1));
            mean_reference_offset2 = rad2deg(circ_mean(visual_offset2(blockLimits{last_high_contrast}(1):blockLimits{last_high_contrast}(2)),[],1));
            
            %% Calculate and plot bump magnitude in time
            
            %1) Get bump magnitude per block
            for block = 1:length(blockLimits)
                new_bump_mag{block} = continuous_data.bump_magnitude(:,blockLimits{block}(1):blockLimits{block}(2));
                adj_rs{block} = continuous_data.adj_rs(:,blockLimits{block}(1):blockLimits{block}(2));
            end
            
            %2) Plot
            figure('Position',[200 200 1600 600]),
            %Plot EPG activity
            subplot(2,1,1)
            imagesc(continuous_data.dff_matrix')
            colormap(flipud(gray))
            hold on
            %add the changes in stim
            for change = 1:length(changeContrast)
                line([changeContrast(change) changeContrast(change)], [0 size(continuous_data.dff_matrix,2)], 'LineWidth', 2, 'color', [0,0.5,0]);
            end
            title('EPG activity in the PB','fontweight','bold','fontsize',12);
            set(gca,'xtick',[]);
            set(gca,'XTickLabel',[]);
            legend('Change in stimulus');
            
            %Plot bump magnitude
            subplot(2,1,2)
            for block = 1:length(blockLimits)
                time = continuous_data.time(blockLimits{block}(1):blockLimits{block}(2));
                if contains(contrasts(block),'Dark')
                    plot(time(adj_rs{block}>=0.5),new_bump_mag{block}(adj_rs{block}>=0.5),'.','color',color_gradient{1})
                    hold on
                    plot(time(adj_rs{block}<0.5),new_bump_mag{block}(adj_rs{block}<0.5),'.r')
                elseif contains(contrasts(block),'Low')
                    plot(time(adj_rs{block}>=0.5),new_bump_mag{block}(adj_rs{block}>=0.5),'.','color',color_gradient{2})
                    hold on
                    plot(time(adj_rs{block}<0.5),new_bump_mag{block}(adj_rs{block}<0.5),'.r')
                else
                    plot(time(adj_rs{block}>=0.5),new_bump_mag{block}(adj_rs{block}>=0.5),'.','color',color_gradient{3})
                    hold on
                    plot(time(adj_rs{block}<0.5),new_bump_mag{block}(adj_rs{block}<0.5),'.r')
                end
            end
            title('Bump magnitude','fontweight','bold','fontsize',12);
            ylabel({'Bump magnitude';'(from von Mises fit)'},'fontweight','bold','fontsize',10);
            xlabel('Time (sec)','fontweight','bold','fontsize',10);
            xlim([0 continuous_data.time(end)]);
            
            %4) Save figure
            saveas(gcf,[path,'continuous_plots\closedLoopBMinTime.png']);
            
            %% Compute and plot mean bump magnitude per block
            
            figure('Position',[200 200 1000 800]),
            subplot(2,1,1)
            mean_new_bump_mag = [];
            for block = 1:length(blockLimits)
                mean_new_bump_mag(block) = mean(new_bump_mag{block}(adj_rs{block}>=0.5)); %select only datapoints with good gof
                plot(block,mean_new_bump_mag(block),'ko','MarkerFaceColor',color_gradient{Intensities(block)},'MarkerSize',8)
                hold on
            end
            xlim([0 num_subplots+1]);
            ylim([min(mean_new_bump_mag)-0.5 max(mean_new_bump_mag)+0.5]);
            xticks(1:num_subplots);
            title('Mean bump magnitude per block');
            ylabel({'Bump magnitude';'(from von Mises fit)'});
            xlabel('Block number');
            %Add custom legend
            h = zeros(3, 1);
            h(1) = plot(NaN,NaN,'ko','MarkerFaceColor',color_gradient{1});
            h(2) = plot(NaN,NaN,'ko','MarkerFaceColor',color_gradient{2});
            h(3) = plot(NaN,NaN,'ko','MarkerFaceColor',color_gradient{3});
            legend(h, 'Darkness','Low contrast','High Contrast');
            
            subplot(2,1,2)
            %Add bump mag data to table
            if size(summary_data,2) == 2
                summary_data = addvars(summary_data, mean_new_bump_mag','NewVariableNames','mean_new_bump_mag');
            end
            %Get mean bump mag by contrast
            mean_new_bump_data = varfun(@mean,summary_data,'InputVariables','mean_new_bump_mag',...
                'GroupingVariables',{'contrast'});
            %Plot
            for contrastLevel = 1:3
                plot(contrastLevel,mean_new_bump_data.mean_mean_new_bump_mag(newOrder(contrastLevel)),'ko','MarkerFaceColor',color_gradient{contrastLevel},'MarkerSize',8)
                hold on
            end
            xlim([0 4]);
            ylim([min(mean_new_bump_data.mean_mean_new_bump_mag)-0.5 max(mean_new_bump_data.mean_mean_new_bump_mag)+0.5]);
            xticks(1:3);
            xticklabels({'Darkness', 'Low contrast', 'High contrast'});
            title('Mean bump magnitude per contrast');
            ylabel({'Bump magnitude';'(from von Mises fit)'});
            
            %save figure
            saveas(gcf,[path,'continuous_plots\closedLoopMeanBMNewMethod.png']);
            
            %% Analyze bump width at half max per block with the new method
            
            %Compute bump width at half maximum (using aux function)
            new_half_max_width = continuous_data.bump_width;
            
            for block = 1:length(blockLimits)
                new_width_half_max{block} = new_half_max_width(blockLimits{block}(1):blockLimits{block}(2));
            end
            
            %Plot
            figure('Position',[200 200 1000 800]),
            subplot(2,1,1)
            for block = 1:length(blockLimits)
                mean_new_half_w(block) = nanmean(new_width_half_max{block}(adj_rs{block}>=0.5));
                plot(block,mean_new_half_w(block),'ko','MarkerFaceColor',color_gradient{Intensities(block)},'MarkerSize',8)
                hold on
            end
            xlim([0 num_subplots+1]);
            xticks(1:num_subplots);
            ylim([min(mean_new_half_w)-0.5 max(mean_new_half_w)+0.5]);
            title('Mean half max width per block');
            ylabel('Half max width (EB wedges)');
            xlabel('Block number');
            %Add custom legend
            h = zeros(3, 1);
            h(1) = plot(NaN,NaN,'ko','MarkerFaceColor',color_gradient{1});
            h(2) = plot(NaN,NaN,'ko','MarkerFaceColor',color_gradient{2});
            h(3) = plot(NaN,NaN,'ko','MarkerFaceColor',color_gradient{3});
            legend(h, 'Darkness','Low contrast','High Contrast');
            
            subplot(2,1,2)
            %Add bump mag data to table
            if size(summary_data,2) == 3
                summary_data = addvars(summary_data, mean_new_half_w','NewVariableNames','mean_new_half_width');
            end
            %Get mean bump mag by contrast
            mean_new_half_width = varfun(@mean,summary_data,'InputVariables','mean_new_half_width',...
                'GroupingVariables',{'contrast'});
            %Plot
            for contrastLevel = 1:3
                plot(contrastLevel,mean_new_half_width.mean_mean_new_half_width(newOrder(contrastLevel)),'ko','MarkerFaceColor',color_gradient{contrastLevel},'MarkerSize',8)
                hold on
            end
            xlim([0 4]);
            ylim([min(mean_new_half_width.mean_mean_new_half_width)-0.5 max(mean_new_half_width.mean_mean_new_half_width)+0.5]);
            xticks(1:3);
            xticklabels({'Darkness', 'Low contrast', 'High contrast'});
            title('Width at half max per contrast');
            ylabel('Full width at half max');
            
            %Save figure
            saveas(gcf,[path,'continuous_plots\closedLoopMeanHWNewMethod.png']);
            
            
            %% Repeat using new method
            
            allNewBumpMag = [];
            all_adj_rs = [];
            for block = 1:length(blockLimits)
                allNewBumpMag = [allNewBumpMag,new_bump_mag{block}];
                all_adj_rs = [all_adj_rs,adj_rs{block}];
            end
            
            %Create vector with the contrast level for each timepoint
            all_contrast_levels = [];
            for block = 1:length(blockLimits)
                contrast_level{block} = repelem(Intensities(block),blockLimits{block}(2)+1-blockLimits{block}(1));
                all_contrast_levels = [all_contrast_levels,contrast_level{block}];
            end
            
            figure('Position',[200 200 1400 600]),
            nbins = 20;
            
            %Forward velocity
            subplot(1,4,1)
            %Define bin limits
            maxBin = max(continuous_data.vel_for_deg_ds); %upper limit
            binWidth = maxBin/nbins;
            forVelBins = [0:binWidth:maxBin];
            
            %Create axes for plot, centering them in the middle of the bin
            forVelAxes = forVelBins-binWidth/2;
            forVelAxes = forVelAxes(2:end);
            
            %Get binned means
            for contrast = 1:3
                for bin = 1:length(forVelBins)-1
                    doubleBin(bin,contrast) = nanmean(allNewBumpMag((continuous_data.vel_for_deg_ds(1:length(allNewBumpMag)) > forVelBins(bin)) & (continuous_data.vel_for_deg_ds(1:length(allNewBumpMag)) < forVelBins(bin+1)) & (all_adj_rs' >= 0.5) & (all_contrast_levels' == contrast)));
                end
                plot(forVelAxes,doubleBin(:,contrast),'-o','color',color_gradient{contrast})
                hold on
            end
            ylabel('Mean bump magnitude'); xlabel('Forward velocity (deg/s)');
            ylim([0 (max(max(doubleBin))+0.5)]);
            %if you have a row with nans, make that the upper limit
            nan_row = sum(isnan(doubleBin(:,:)),2);
            nan_rows = find(nan_row==3);
            if ~isempty(nan_rows)
                xlim([0 forVelAxes(nan_rows(1)-1)+binWidth/2]);
            else
                xlim([0 forVelAxes(end)+binWidth/2]);
            end
            legend('Darkness','Low contrast','High contrast');
            
            
            %Side speed
            subplot(1,4,2)
            
            %Define bin limits
            maxBin = max(abs(continuous_data.vel_side_deg_ds));
            binWidth = maxBin/nbins;
            sideSpeedBins = [0:binWidth:maxBin];
            
            %Create axes for plot
            sideSpeedAxes = sideSpeedBins-binWidth/2;
            sideSpeedAxes = sideSpeedAxes(2:end);
            
            %Get binned means
            for contrast = 1:3
                for bin = 1:length(sideSpeedBins)-1
                    doubleBin(bin,contrast) = nanmean(allNewBumpMag((abs(continuous_data.vel_side_deg_ds(1:length(allNewBumpMag))) > sideSpeedBins(bin)) & (abs(continuous_data.vel_side_deg_ds(1:length(allNewBumpMag))) < sideSpeedBins(bin+1)) & (all_adj_rs' >= 0.5) & (all_contrast_levels' == contrast)));
                end
                plot(sideSpeedAxes,doubleBin(:,contrast),'-o','color',color_gradient{contrast})
                hold on
            end
            ylabel('Mean bump magnitude'); xlabel('Side speed (deg/s)');
            ylim([0 (max(max(doubleBin))+0.5)]);
            %if you have a row with nans, make that the upper limit
            nan_row = sum(isnan(doubleBin(:,:)),2);
            nan_rows = find(nan_row==3);
            if ~isempty(nan_rows)
                xlim([0 sideSpeedAxes(nan_rows(1)-1)+binWidth/2]);
            else
                xlim([0 sideSpeedAxes(end)+binWidth/2]);
            end
            legend('Darkness','Low contrast','High contrast');
            
            
            %Angular speed
            subplot(1,4,3)
            
            %Define bin limits
            maxBin = max(abs(continuous_data.vel_yaw_ds));
            binWidth = maxBin/nbins;
            angSpeedBins = [0:binWidth:maxBin];
            
            %Create axes for plot
            angSpeedAxes = angSpeedBins-binWidth/2;
            angSpeedAxes = angSpeedAxes(2:end);
            
            %Get binned means
            for contrast = 1:3
                for bin = 1:length(angSpeedBins)-1
                    doubleBin(bin,contrast) = nanmean(allNewBumpMag((abs(continuous_data.vel_yaw_ds(1:length(allNewBumpMag))) > angSpeedBins(bin)) & (abs(continuous_data.vel_yaw_ds(1:length(allNewBumpMag))) < angSpeedBins(bin+1)) & (all_adj_rs >= 0.5) & (all_contrast_levels == contrast)));
                end
                plot(angSpeedAxes,doubleBin(:,contrast),'-o','color',color_gradient{contrast})
                hold on
            end
            ylabel('Mean bump magnitude'); xlabel('Angular speed (deg/s)');
            ylim([0 (max(max(doubleBin))+0.5)]);
            %if you have a row with nans, make that the upper limit
            nan_row = sum(isnan(doubleBin(:,:)),2);
            nan_rows = find(nan_row==3);
            if ~isempty(nan_rows)
                xlim([0 angSpeedAxes(nan_rows(1)-1)+binWidth/2]);
            else
                xlim([0 angSpeedAxes(end)+binWidth/2]);
            end
            legend('Darkness','Low contrast','High contrast');
            
            
            
            %Total movement
            subplot(1,4,4)
            
            %Define bin limits
            maxBin = max(continuous_data.total_mvt_ds);
            binWidth = maxBin/nbins;
            mvtBins = [0:binWidth:maxBin];
            
            %Create axes for plot
            mvtAxes = mvtBins-binWidth/2;
            mvtAxes = mvtAxes(2:end);
            
            %Get binned means
            for contrast = 1:3
                for bin = 1:length(mvtBins)-1
                    doubleBin(bin,contrast) = mean(allNewBumpMag((continuous_data.total_mvt_ds(1:length(allNewBumpMag)) > mvtBins(bin)) & (continuous_data.total_mvt_ds(1:length(allNewBumpMag)) < mvtBins(bin+1)) & (all_adj_rs >= 0.5) & (all_contrast_levels == contrast)));
                end
                plot(mvtAxes,doubleBin(:,contrast),'-o','color',color_gradient{contrast})
                hold on
            end
            ylabel('Mean bump magnitude'); xlabel('Total movement (deg/s)');
            %if you have a row with nans, make that the upper limit
            nan_row = sum(isnan(doubleBin(:,:)),2);
            nan_rows = find(nan_row==3);
            if ~isempty(nan_rows)
                xlim([0 mvtAxes(nan_rows(1)-1)+binWidth/2]);
            else
                xlim([0 mvtAxes(end)+binWidth/2]);
            end
            ylim([0 (max(max(doubleBin))+0.5)]);
            legend('Darkness','Low contrast','High contrast');
            
            %save figure
            saveas(gcf,[path,'continuous_plots\closedLoopBMVelPerContrastBinned.png']);
            
            %% Model bump magnitude as a function of contrast level and total movement
            
            %fill missing data in the movement variables
            continuous_data.vel_for_deg_ds = fillmissing(continuous_data.vel_for_deg_ds,'linear');
            continuous_data.vel_side_deg_ds = fillmissing(continuous_data.vel_side_deg_ds,'linear');
            continuous_data.vel_yaw_ds = fillmissing(continuous_data.vel_yaw_ds,'linear');
            continuous_data.total_mvt_ds = fillmissing(continuous_data.total_mvt_ds,'linear');
            %zscore the movement data
            zscored_for_vel = zscore(continuous_data.vel_for_deg_ds);
            zscored_side_speed = zscore(abs(continuous_data.vel_side_deg_ds));
            zscored_yaw_speed = zscore(abs(continuous_data.vel_yaw_ds));
            zscored_total_mvt = zscore(continuous_data.total_mvt_ds);
            
            %Create table with the model's variables
            modelTable = table(all_contrast_levels',continuous_data.vel_for_deg_ds(1:length(allNewBumpMag)),zscored_for_vel(1:length(allNewBumpMag)),abs(continuous_data.vel_side_deg_ds(1:length(allNewBumpMag))),zscored_side_speed(1:length(allNewBumpMag)),abs(continuous_data.vel_yaw_ds(1:length(allNewBumpMag)))',zscored_yaw_speed(1:length(allNewBumpMag))',continuous_data.total_mvt_ds(1:length(allNewBumpMag))',zscored_total_mvt(1:length(allNewBumpMag))',continuous_data.time(1:length(allNewBumpMag)),allNewBumpMag',all_adj_rs','VariableNames',{'ContrastLevel','ForVelocity','ZscoredForVel','SideSpeed','ZscoredSideSpeed','YawSpeed','ZscoredYawSpeed','TotalMovement','ZscoredTotalMovement','Time','NewBumpMagnitude','AdjRSquare'});
            
            %% Relationship between bump width at half max and velocity
            
            allNewBumpWidth = [];
            for block = 1:length(blockLimits)
                allNewBumpWidth = [allNewBumpWidth,new_width_half_max{block}];
            end
            
            %Plot
            figure('Position',[200 200 1400 600]),
            
            %Forward velocity
            subplot(1,4,1)
            
            %Define bin limits
            maxBin = max(continuous_data.vel_for_deg_ds); %upper limit
            binWidth = maxBin/nbins;
            forVelBins = [0:binWidth:maxBin];
            
            %Create axes for plot, centering them in the middle of the bin
            forVelAxes = forVelBins-binWidth/2;
            forVelAxes = forVelAxes(2:end);
            
            %Get binned means
            for contrast = 1:3
                for bin = 1:length(forVelBins)-1
                    doubleBin(bin,contrast) = nanmean(allNewBumpWidth((continuous_data.vel_for_deg_ds(1:length(allNewBumpWidth)) > forVelBins(bin)) & (continuous_data.vel_for_deg_ds(1:length(allNewBumpWidth)) < forVelBins(bin+1)) & (all_adj_rs' >= 0.5) & (all_contrast_levels' == contrast)));
                end
                plot(forVelAxes,doubleBin(:,contrast),'-o','color',color_gradient{contrast})
                hold on
            end
            ylabel('Mean bump width'); xlabel('Forward velocity (deg/s)');
            ylim([0 (max(max(doubleBin))+0.5)]);
            %if you have a row with nans, make that the upper limit
            nan_row = sum(isnan(doubleBin(:,:)),2);
            nan_rows = find(nan_row==3);
            if ~isempty(nan_rows)
                xlim([0 forVelAxes(nan_rows(1)-1)+binWidth/2]);
            else
                xlim([0 forVelAxes(end)+binWidth/2]);
            end
            legend('Darkness','Low contrast','High contrast');
            
            
            
            %Side speed
            subplot(1,4,2)
            
            %Define bin limits
            maxBin = max(abs(continuous_data.vel_side_deg_ds));
            binWidth = maxBin/nbins;
            sideSpeedBins = [0:binWidth:maxBin];
            
            %Create axes for plot
            sideSpeedAxes = sideSpeedBins-binWidth/2;
            sideSpeedAxes = sideSpeedAxes(2:end);
            
            %Get binned means
            for contrast = 1:3
                for bin = 1:length(sideSpeedBins)-1
                    doubleBin(bin,contrast) = nanmean(allNewBumpWidth((abs(continuous_data.vel_side_deg_ds(1:length(allNewBumpWidth))) > sideSpeedBins(bin)) & (abs(continuous_data.vel_side_deg_ds(1:length(allNewBumpWidth))) < sideSpeedBins(bin+1)) & (all_adj_rs' >= 0.5) & (all_contrast_levels' == contrast)));
                end
                plot(sideSpeedAxes,doubleBin(:,contrast),'-o','color',color_gradient{contrast})
                hold on
            end
            ylabel('Mean bump width'); xlabel('Side speed (deg/s)');
            ylim([0 (max(max(doubleBin))+0.5)]);
            %if you have a row with nans, make that the upper limit
            nan_row = sum(isnan(doubleBin(:,:)),2);
            nan_rows = find(nan_row==3);
            if ~isempty(nan_rows)
                xlim([0 sideSpeedAxes(nan_rows(1)-1)+binWidth/2]);
            else
                xlim([0 sideSpeedAxes(end)+binWidth/2]);
            end
            legend('Darkness','Low contrast','High contrast');
            
            
            %Angular speed
            subplot(1,4,3)
            
            %Define bin limits
            maxBin = max(abs(continuous_data.vel_yaw_ds));
            binWidth = maxBin/nbins;
            angSpeedBins = [0:binWidth:maxBin];
            
            %Create axes for plot
            angSpeedAxes = angSpeedBins-binWidth/2;
            angSpeedAxes = angSpeedAxes(2:end);
            
            %Get binned means
            for contrast = 1:3
                for bin = 1:length(angSpeedBins)-1
                    doubleBin(bin,contrast) = nanmean(allNewBumpWidth((abs(continuous_data.vel_yaw_ds(1:length(allNewBumpWidth))) > angSpeedBins(bin)) & (abs(continuous_data.vel_yaw_ds(1:length(allNewBumpWidth))) < angSpeedBins(bin+1)) & (all_adj_rs >= 0.5) & (all_contrast_levels == contrast)));
                end
                plot(angSpeedAxes,doubleBin(:,contrast),'-o','color',color_gradient{contrast})
                hold on
            end
            ylabel('Mean bump width'); xlabel('Angular speed (deg/s)');
            ylim([0 (max(max(doubleBin))+0.5)]);
            %if you have a row with nans, make that the upper limit
            nan_row = sum(isnan(doubleBin(:,:)),2);
            nan_rows = find(nan_row==3);
            if ~isempty(nan_rows)
                xlim([0 angSpeedAxes(nan_rows(1)-1)+binWidth/2]);
            else
                xlim([0 angSpeedAxes(end)+binWidth/2]);
            end
            legend('Darkness','Low contrast','High contrast');
            
            
            
            %Total movement
            subplot(1,4,4)
            
            %Define bin limits
            maxBin = max(continuous_data.total_mvt_ds);
            binWidth = maxBin/nbins;
            mvtBins = [0:binWidth:maxBin];
            
            %Create axes for plot
            mvtAxes = mvtBins-binWidth/2;
            mvtAxes = mvtAxes(2:end);
            
            %Get binned means
            for contrast = 1:3
                for bin = 1:length(mvtBins)-1
                    doubleBin(bin,contrast) = mean(allNewBumpWidth((continuous_data.total_mvt_ds(1:length(allNewBumpWidth)) > mvtBins(bin)) & (continuous_data.total_mvt_ds(1:length(allNewBumpWidth)) < mvtBins(bin+1)) & (all_adj_rs >= 0.5) & (all_contrast_levels == contrast)));
                end
                plot(mvtAxes,doubleBin(:,contrast),'-o','color',color_gradient{contrast})
                hold on
            end
            ylabel('Mean bump width'); xlabel('Total movement (deg/s)');
            %if you have a row with nans, make that the upper limit
            nan_row = sum(isnan(doubleBin(:,:)),2);
            nan_rows = find(nan_row==3);
            if ~isempty(nan_rows)
                xlim([0 mvtAxes(nan_rows(1)-1)+binWidth/2]);
            else
                xlim([0 mvtAxes(end)+binWidth/2]);
            end
            ylim([0 (max(max(doubleBin))+0.5)]);
            legend('Darkness','Low contrast','High contrast');
            
            %save figure
            saveas(gcf,[path,'continuous_plots\closedLoopHWNewMethodvsVelPerContrastBinned.png']);
            
            %% Model bump width at half max as a function of contrastLevel and total movement
            
            %Create table with the model's variables
            modelTable = addvars(modelTable,allNewBumpWidth','NewVariableNames','NewBumpWidth');

            %% Heading variation per stimulus
            
            figure('Position',[200 200 1000 800]),
            subplot(2,1,1)
            for block = 1:length(blockLimits)
                [~, heading_var(block)] = circ_std(deg2rad(heading(blockLimits{block}(1):blockLimits{block}(2))),[],[],1);
                plot(block,heading_var(block),'ko','MarkerFaceColor',color_gradient{Intensities(block)},'MarkerSize',8)
                hold on
            end
            xlim([0 num_subplots+1]);
            ylim([min(heading_var) - 0.5 max(heading_var) + 0.5]);
            xticks(1:num_subplots);
            title('Heading variation per block');
            ylabel({'Circular standard deviation','of the heading'});
            xlabel('Block number');
            %Add custom legend
            h = zeros(3, 1);
            h(1) = plot(NaN,NaN,'ko','MarkerFaceColor',color_gradient{1});
            h(2) = plot(NaN,NaN,'ko','MarkerFaceColor',color_gradient{2});
            h(3) = plot(NaN,NaN,'ko','MarkerFaceColor',color_gradient{3});
            legend(h, 'Darkness','Low contrast','High Contrast');
            
            subplot(2,1,2)
            %Create table with contrast level and offset variation
            if size(summary_data,2) == 4
                summary_data = addvars(summary_data, heading_var','NewVariableNames','heading_var');
            end
            %Get mean heading var by contrast
            mean_heading_data = varfun(@mean,summary_data,'InputVariables','heading_var',...
                'GroupingVariables',{'contrast'});
            %Plot
            for contrastLevel = 1:3
                plot(contrastLevel,mean_heading_data.mean_heading_var(newOrder(contrastLevel)),'ko','MarkerFaceColor',color_gradient{contrastLevel},'MarkerSize',8)
                hold on
            end
            xlim([0 4]);
            ylim([min(mean_heading_data.mean_heading_var)-0.5 max(mean_heading_data.mean_heading_var)+0.5]);
            xticks(1:3);
            xticklabels({'Darkness', 'Low contrast', 'High contrast'});
            title('Mean heading variation per contrast');
            ylabel({'Circular standard deviation','of the heading'});
            
            %save figure
            saveas(gcf,[path,'continuous_plots\closedLoopHeadingVariation.png']);
            
            
            %% Save useful data
            
            save([path,'continuous_summary_data.mat'],'summary_data','mean_reference_offset','mean_reference_offset2','modelTable')
            
                       
            clearvars -except data_dirs fly fly_files sid
            close all
            
        else
        end
    end
    
end
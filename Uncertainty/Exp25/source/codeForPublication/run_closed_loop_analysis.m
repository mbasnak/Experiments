
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
            
            
            %% If this is fly 2, remove all datapoints beyond 8000, since the panels malfunctionned there
            
            if (contains(path, '20201020_60D05_7f') & ~contains(path,'fly2'))
                data.dff_pva = data.dff_pva(1:8000);
                data.heading_deg = data.heading_deg(1:8000);
                data.phase = data.phase(1:8000);
                data.offset = data.offset(1:8000);
                data.time = data.time(1:8000);
                data.dff_matrix = data.dff_matrix(:,1:8000);
                data.mean_dff_EB = data.mean_dff_EB(:,1:8000);
                data.fr_y_ds = data.fr_y_ds(1:8000);
                data.vel_for_deg_ds = data.vel_for_deg_ds(1:8000);
                data.vel_side_deg_ds = data.vel_side_deg_ds(1:8000);
                data.vel_yaw_ds = data.vel_yaw_ds(1:8000);
                data.total_mvt_ds = data.total_mvt_ds(1:8000);
                data.panel_angle = data.panel_angle(1:8000);
            end
            %% Identify changes in stimulus
            
            %Plot the signal from the y dimension of the panels
            figure,
            subplot 121,
            plot(data.fr_y_ds)
            title('Panels y signal');
            xlabel('Time (downsampled frames)');
            ylabel('Voltage');
            xlim([0 length(data.fr_y_ds)]);
            
            subplot 122,
            change_y_panels = abs(diff(data.fr_y_ds));
            plot(change_y_panels)
            %Obtain the frames where the stim changes using the derivative of the
            %signal from the y dimension of the panels
            changeContrast = find(abs(diff(data.fr_y_ds))>1);
            hold on
            %Add them to the plot
            plot(changeContrast,change_y_panels(changeContrast),'ro')
            title('Change in y panels');
            xlim([0 length(data.fr_y_ds)]);
            xlabel('Time (downsampled frames)');
            
            
            %% Identify contrast order
            
            %Recover the position function we used to change the pattern
            pos_function = data.run_obj.function_number;
            
            %Define the order of intensities presented according to the function used
            if pos_function == 195
                Intensities = [1,2,1,3,2,3];
            elseif pos_function == 196
                Intensities = [2,1,3,1,2,3];
            else
                Intensities = [3,1,2,1,2,3];
            end
            
            if (contains(path, '20201020_60D05_7f') & ~contains(path,'fly2'))
                Intensities = Intensities(1:5);
            end
            %% Plot the activity heatmap, phase and heading positions, and offset in time
            
            %Define the number of columns in the subplots depending on the fly (the one
            %where fictrac malfunctioned will have one less column)
            if (contains(path, '20201020_60D05_7f') & ~contains(path,'fly2'))
                num_subplots = 5;
            else
                num_subplots = 6;
            end
            
            %Plot
            % Plot the heatmap of EPG activity
            figure('Position',[0 0 1800 800]),
            subplot(4,num_subplots,[1 num_subplots])
            %I'm flipping the dff matrix for it to make sense along with the fly's
            %heading
            imagesc(flip(data.dff_matrix))
            colormap(flipud(gray))
            hold on
            %add the changes in stim
            for change = 1:length(changeContrast)
                line([changeContrast(change) changeContrast(change)], [0 17], 'LineWidth', 3, 'color', [0, 0.5, 0]);
            end
            yticks(1:2:16);
            yticklabels({'8L','6L','4L','2L','1R','3R','5R','7R'});
            ylabel('PB glomerulus','fontweight','bold','fontsize',10);
            title('EPG activity in the PB','fontweight','bold','fontsize',12);
            set(gca,'XTickLabel',[]);
            legend('Change in stimulus');
            
            % Plot the heading and the EPG phase
            subplot(4,num_subplots,[num_subplots+1 num_subplots*2])
            %Get heading to plot
            heading = wrapTo180(data.heading_deg);
            %Remove wrapped lines to plot
            [x_out_heading,heading_to_plot] = removeWrappedLines(data.time,heading);
            plot(x_out_heading,heading_to_plot,'LineWidth',1.5,'color',[0.2 0.6 0.7])
            hold on
            %Get EPG phase to plot
            %I'm now going to negate the phase, since I'm plotting heading instead of
            %bar position, so the bump moves in the other direction
            phase = wrapTo180(rad2deg(-data.dff_pva));
            [x_out_phase,phase_to_plot] = removeWrappedLines(data.time,phase);
            plot(x_out_phase,phase_to_plot,'color',[0.9 0.3 0.4],'LineWidth',1.5)
            %add the changes in stim
            for change = 1:length(changeContrast)
                line([data.time(changeContrast(change)) data.time(changeContrast(change))], [-180 180], 'LineWidth', 3, 'color', [0, 0.5, 0]);
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
            offset = wrapTo180(data.offset);
            [x_out_offset,offset_to_plot] = removeWrappedLines(data.time,offset);
            plot(x_out_offset,offset_to_plot,'LineWidth',1.5,'color','k')
            %Add the changes in stim
            for change = 1:length(changeContrast)
                line([data.time(changeContrast(change)) data.time(changeContrast(change))], [-180 180], 'LineWidth', 3, 'color', [0, 0.5, 0]);
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
            saveas(gcf,[path,'plots\closedLoopHeatmapAndOffset.png']);
            
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
            legend(h, 'Darkness','Low contrast','High Contrast');
            
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
            saveas(gcf,[path,'plots\closedLoopOffsetVariation.png']);
            
            %% Get offset value from last bout of high contrast: this will be our 'reference offset' to use during the open loop analysis
            
            %Find limits of last high contrast bout
            high_contrast_bouts = find(Intensities(1:length(blockLimits))==3);
            last_high_contrast = high_contrast_bouts(end);
            
            %This will be the 'visual' offset, so we need to compute it as the circular
            %distance between the stimulus and the phase
            visual_offset = circ_dist(data.dff_pva,deg2rad(data.panel_angle));
            
            %Compute the circular mean of this visual offset in the last high contrast
            %bout
            mean_reference_offset = rad2deg(circ_mean(visual_offset(blockLimits{last_high_contrast}(1):blockLimits{last_high_contrast}(2)),[],1));
            
            %% Calculate and plot bump magnitude in time
            
            %Get bump magnitude per block
            for block = 1:length(blockLimits)
                bump_mag{block} = data.bump_magnitude(:,blockLimits{block}(1):blockLimits{block}(2));
            end
            
            figure('Position',[200 200 1600 600]),
            %Plot EPG activity
            subplot(2,1,1)
            imagesc(flip(data.dff_matrix))
            colormap(gray)
            hold on
            %add the changes in stim
            for change = 1:length(changeContrast)
                line([changeContrast(change) changeContrast(change)], [0 17], 'LineWidth', 2, 'color', [0,0.5,0]);
            end
            yticks(1:2:16);
            yticklabels({'8L','6L','4L','2L','1R','3R','5R','7R'});
            ylabel('PB glomerulus','fontweight','bold','fontsize',10);
            title('EPG activity in the PB','fontweight','bold','fontsize',12);
            set(gca,'xtick',[]);
            set(gca,'XTickLabel',[]);
            legend('Change in stimulus');
            
            %Plot bump magnitude
            subplot(2,1,2)
            for block = 1:length(blockLimits)
                if contains(contrasts(block),'Dark')
                    plot(data.time(blockLimits{block}(1):blockLimits{block}(2)),bump_mag{block},'color',color_gradient{1})
                elseif contains(contrasts(block),'Low')
                    plot(data.time(blockLimits{block}(1):blockLimits{block}(2)),bump_mag{block},'color',color_gradient{2})
                else
                    plot(data.time(blockLimits{block}(1):blockLimits{block}(2)),bump_mag{block},'color',color_gradient{3})
                end
                hold on
            end
            title('Bump magnitude','fontweight','bold','fontsize',12);
            ylabel({'Bump magnitude';'(amplitude of Fourier component)'},'fontweight','bold','fontsize',10);
            xlabel('Time (sec)','fontweight','bold','fontsize',10);
            xlim([0 data.time(end)]);
            
            %save figure
            saveas(gcf,[path,'plots\closedLoopBMinTime.png']);
            
            %% Compute and plot mean bump magnitude per block
            
            figure('Position',[200 200 1000 800]),
            subplot(2,1,1)
            for block = 1:length(blockLimits)
                mean_bump_mag(block) = mean(bump_mag{block});
                plot(block,mean_bump_mag(block),'ko','MarkerFaceColor',color_gradient{Intensities(block)},'MarkerSize',8)
                hold on
            end
            xlim([0 num_subplots+1]);
            ylim([min(mean_bump_mag)-0.5 max(mean_bump_mag)+0.5]);
            xticks(1:num_subplots);
            title('Mean bump magnitude per block');
            ylabel({'Bump magnitude';'(amplitude of Fourier component)'});
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
                summary_data = addvars(summary_data, mean_bump_mag','NewVariableNames','mean_bump_mag');
            end
            %Get mean bump mag by contrast
            mean_bump_data = varfun(@mean,summary_data,'InputVariables','mean_bump_mag',...
                'GroupingVariables',{'contrast'});
            %Plot
            for contrastLevel = 1:3
                plot(contrastLevel,mean_bump_data.mean_mean_bump_mag(newOrder(contrastLevel)),'ko','MarkerFaceColor',color_gradient{contrastLevel},'MarkerSize',8)
                hold on
            end
            xlim([0 4]);
            ylim([min(mean_bump_data.mean_mean_bump_mag)-0.5 max(mean_bump_data.mean_mean_bump_mag)+0.5]);
            xticks(1:3);
            xticklabels({'Darkness', 'Low contrast', 'High contrast'});
            title('Mean bump magnitude per contrast');
            ylabel({'Bump magnitude';'(amplitude of Fourier component)'});
            
            %save figure
            saveas(gcf,[path,'plots\closedLoopMeanBM.png']);
            
            %% Compute bump width at half maximum
            
            for timepoint = 1:size(data.mean_dff_EB,2)
                half_max_width(timepoint) = compute_bump_width(data.mean_dff_EB,timepoint);
            end
            
            %% Analyze bump width at half max per block
            
            for block = 1:length(blockLimits)
                width_half_max{block} = half_max_width(blockLimits{block}(1):blockLimits{block}(2));
            end
            
            %plot
            figure('Position',[200 200 1000 800]),
            subplot(2,1,1)
            for block = 1:length(blockLimits)
                mean_half_w(block) = nanmean(width_half_max{block});
                plot(block,mean_half_w(block),'ko','MarkerFaceColor',color_gradient{Intensities(block)},'MarkerSize',8)
                hold on
            end
            xlim([0 num_subplots+1]);
            xticks(1:num_subplots);
            ylim([min(mean_half_w)-0.5 max(mean_half_w)+0.5]);
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
                summary_data = addvars(summary_data, mean_half_w','NewVariableNames','mean_half_width');
            end
            %Get mean bump mag by contrast
            mean_half_width = varfun(@mean,summary_data,'InputVariables','mean_half_width',...
                'GroupingVariables',{'contrast'});
            %Plot
            for contrastLevel = 1:3
                plot(contrastLevel,mean_half_width.mean_mean_half_width(newOrder(contrastLevel)),'ko','MarkerFaceColor',color_gradient{contrastLevel},'MarkerSize',8)
                hold on
            end
            xlim([0 4]);
            ylim([min(mean_half_width.mean_mean_half_width)-0.5 max(mean_half_width.mean_mean_half_width)+0.5]);
            xticks(1:3);
            xticklabels({'Darkness', 'Low contrast', 'High contrast'});
            title('Width at half max per contrast');
            ylabel('Full width at half max');
            
            %save figure
            saveas(gcf,[path,'plots\closedLoopMeanHW.png']);
            
            %% Plot fly's velocity in all 3 axes  with all vel in deg/s
            
            figure('Position',[100, 100, 1600, 1000]),
            
            %Forward velocity
            subplot(4,5,[1 4])
            plot(data.time,data.vel_for_deg_ds)
            hold on
            yline(0,'--');
            title('Forward velocity');
            ylabel('Forward velocity (deg/s)');
            xlim([0 data.time(end)]);
            
            subplot(4,5,5)
            histogram(data.vel_for_deg_ds)
            title('Velocity distributions');
            xlabel('Forward velocity (deg/s)');
            ylabel('Counts');
            
            %Side velocity
            subplot(4,5,[6 9])
            plot(data.time,data.vel_side_deg_ds,'color',[0.8 0.2 0.6])
            title('Side velocity');
            ylabel('Side velocity (deg/s)');
            hold on
            yline(0,'--');
            xlim([0 data.time(end)]);
            
            subplot(4,5,10)
            histogram(data.vel_side_deg_ds,'FaceColor',[0.8 0.2 0.6])
            xlabel('Side velocity (deg/s)');
            ylabel('Counts');
            
            %Angular velocity
            subplot(4,5,[11 14])
            plot(data.time,data.vel_yaw_ds,'color',[0.4 0.2 0.8])
            title('Angular velocity');
            ylabel('Angular velocity (deg/s)');
            hold on
            yline(0,'--');
            xlim([0 data.time(end)]);
            
            subplot(4,5,15)
            histogram(data.vel_yaw_ds,'FaceColor',[0.4 0.2 0.8])
            xlabel('Angular velocity (deg/s)');
            ylabel('Counts');
            
            %Total movement
            subplot(4,5,[16 19])
            plot(data.time,data.total_mvt_ds,'k')
            title('Total movement');
            ylabel('Total movement (deg/s)');
            xlabel('Time (sec)');
            xlim([0 data.time(end)]);
            
            subplot(4,5,20)
            histogram(data.total_mvt_ds,'FaceColor','k')
            xlabel('Total movement (deg/s)');
            ylabel('Counts');
            
            %save figure
            saveas(gcf,[path,'plots\closedLoopVelDistributions.png']);
            
            %% Compute and plot mean total movement per block
            
            figure('Position',[200 200 1400 800]),
            %Plot total movement in time per block, colored appropriately
            subplot(2,5,[1 3])
            for block = 1:length(blockLimits)
                if contains(contrasts(block),'Dark')
                    plot(data.time(blockLimits{block}(1):blockLimits{block}(2)),data.total_mvt_ds(1,blockLimits{block}(1):blockLimits{block}(2)),'color',color_gradient{1})
                elseif contains(contrasts(block),'Low')
                    plot(data.time(blockLimits{block}(1):blockLimits{block}(2)),data.total_mvt_ds(1,blockLimits{block}(1):blockLimits{block}(2)),'color',color_gradient{2})
                else
                    plot(data.time(blockLimits{block}(1):blockLimits{block}(2)),data.total_mvt_ds(1,blockLimits{block}(1):blockLimits{block}(2)),'color',color_gradient{3})
                end
                hold on
            end
            title('Total movement in time');
            ylabel('Total movement (deg/s)');
            xlabel('Time (sec)');
            xlim([0 data.time(end)]);
            ylim([0 max(data.total_mvt_ds)+10]);
            
            %Plot median total movement per block
            subplot(2,5,[6 8])
            for block = 1:length(blockLimits)
                mean_total_mvt(block) = nanmean(data.total_mvt_ds(1,blockLimits{block}(1):blockLimits{block}(2)));
                plot(block,mean_total_mvt(block),'ko','MarkerFaceColor',color_gradient{Intensities(block)},'MarkerSize',8)
                hold on
            end
            xlim([0 num_subplots+1]);
            xticks(1:num_subplots);
            title('Mean total movement per block');
            ylabel('Total movement (deg/s)');
            ylim([0 max(mean_total_mvt)+10]);
            xlabel('Block #');
            %Add custom legend
            h = zeros(3, 1);
            h(1) = plot(NaN,NaN,'ko','MarkerFaceColor',color_gradient{1});
            h(2) = plot(NaN,NaN,'ko','MarkerFaceColor',color_gradient{2});
            h(3) = plot(NaN,NaN,'ko','MarkerFaceColor',color_gradient{3});
            legend(h, 'Darkness','Low contrast','High Contrast');
            
            %Plot median total movement per block type
            subplot(2,5,[4,5,9,10])
            for contrast = 1:3
                mean_mean_total_mvt(contrast) = mean(mean_total_mvt(Intensities==contrast));
            end
            plot(mean_mean_total_mvt,'-ko','LineWidth',2)
            ylabel('Total movement (deg/s)');
            title('Mean total movement');
            ylim([0 max(mean_mean_total_mvt+10)]);
            xlim([0 4]);
            xticks(1:3);
            xticklabels({'Darkness','Low contrast','High contrast'});
            
            %save figure
            saveas(gcf,[path,'plots\closedLoopTotalMvtPerBlock.png']);
            
            %% Repeat for forward velocity
            
            figure('Position',[200 200 1400 800]),
            %Plot fwd vel in time per block, colored appropriately
            subplot(2,5,[1 3])
            for block = 1:length(blockLimits)
                if contains(contrasts(block),'Dark')
                    plot(data.time(blockLimits{block}(1):blockLimits{block}(2)),data.vel_for_ds(blockLimits{block}(1):blockLimits{block}(2)),'color',color_gradient{1})
                elseif contains(contrasts(block),'Low')
                    plot(data.time(blockLimits{block}(1):blockLimits{block}(2)),data.vel_for_ds(blockLimits{block}(1):blockLimits{block}(2)),'color',color_gradient{2})
                else
                    plot(data.time(blockLimits{block}(1):blockLimits{block}(2)),data.vel_for_ds(blockLimits{block}(1):blockLimits{block}(2)),'color',color_gradient{3})
                end
                hold on
            end
            title('Forward velocity in time');
            ylabel('Foward velocity (mm/s)');
            xlabel('Time (sec)');
            xlim([0 data.time(end)]);
            ylim([0 max(data.vel_for_ds)+1]);
            
            %Plot median total movement per block
            subplot(2,5,[6 8])
            for block = 1:length(blockLimits)
                mean_fwd_vel(block) = nanmean(data.vel_for_ds(blockLimits{block}(1):blockLimits{block}(2)));
                plot(block,mean_fwd_vel(block),'ko','MarkerFaceColor',color_gradient{Intensities(block)},'MarkerSize',8)
                hold on
            end
            xlim([0 num_subplots+1]);
            xticks(1:num_subplots);
            title('Mean forward velocity per block');
            ylabel('Forward velocity (mm/s)');
            ylim([0 max(mean_fwd_vel)+1]);
            xlabel('Block #');
            %Add custom legend
            h = zeros(3, 1);
            h(1) = plot(NaN,NaN,'ko','MarkerFaceColor',color_gradient{1});
            h(2) = plot(NaN,NaN,'ko','MarkerFaceColor',color_gradient{2});
            h(3) = plot(NaN,NaN,'ko','MarkerFaceColor',color_gradient{3});
            legend(h, 'Darkness','Low contrast','High Contrast');
            
            %Plot mean total movement per block type
            subplot(2,5,[4,5,9,10])
            for contrast = 1:3
                mean_mean_fwd_vel(contrast) = mean(mean_fwd_vel(Intensities==contrast));
            end
            plot(mean_mean_fwd_vel,'-ko','LineWidth',2)
            ylabel('Forward velocity (mm/s)');
            title('Mean forward velocity');
            ylim([0 max(mean_mean_fwd_vel)+1]);
            xlim([0 4]);
            xticks(1:3);
            xticklabels({'Darkness','Low contrast','High contrast'});
            
            %save figure
            saveas(gcf,[path,'plots\closedLoopFwdVelPerBlock.png']);
            
            %% Repeat for angular speed
            
            figure('Position',[200 200 1400 800]),
            %Plot ang vel in time per block, colored appropriately
            subplot(2,5,[1 3])
            for block = 1:length(blockLimits)
                if contains(contrasts(block),'Dark')
                    plot(data.time(blockLimits{block}(1):blockLimits{block}(2)),abs(data.vel_yaw_ds(blockLimits{block}(1):blockLimits{block}(2))),'color',color_gradient{1})
                elseif contains(contrasts(block),'Low')
                    plot(data.time(blockLimits{block}(1):blockLimits{block}(2)),abs(data.vel_yaw_ds(blockLimits{block}(1):blockLimits{block}(2))),'color',color_gradient{2})
                else
                    plot(data.time(blockLimits{block}(1):blockLimits{block}(2)),abs(data.vel_yaw_ds(blockLimits{block}(1):blockLimits{block}(2))),'color',color_gradient{3})
                end
                hold on
            end
            title('Angular speed in time');
            ylabel('Angular speed (deg/s)');
            xlabel('Time (sec)');
            xlim([0 data.time(end)]);
            ylim([0 max(abs(data.vel_yaw_ds))+10]);
            
            %Plot median total movement per block
            subplot(2,5,[6 8])
            for block = 1:length(blockLimits)
                mean_ang_speed(block) = nanmean(abs(data.vel_yaw_ds(blockLimits{block}(1):blockLimits{block}(2))));
                plot(block,mean_ang_speed(block),'ko','MarkerFaceColor',color_gradient{Intensities(block)},'MarkerSize',8)
                hold on
            end
            xlim([0 num_subplots+1]);
            xticks(1:num_subplots);
            title('Mean angular speed per block');
            ylabel('Angular speed (deg/s)');
            ylim([0 max(mean_ang_speed)+10]);
            xlabel('Block #');
            %Add custom legend
            h = zeros(3, 1);
            h(1) = plot(NaN,NaN,'ko','MarkerFaceColor',color_gradient{1});
            h(2) = plot(NaN,NaN,'ko','MarkerFaceColor',color_gradient{2});
            h(3) = plot(NaN,NaN,'ko','MarkerFaceColor',color_gradient{3});
            legend(h, 'Darkness','Low contrast','High Contrast');
            
            %Plot median total movement per block type
            subplot(2,5,[4,5,9,10])
            for contrast = 1:3
                mean_mean_ang_speed(contrast) = mean(mean_ang_speed(Intensities==contrast));
            end
            plot(mean_mean_ang_speed,'-ko','LineWidth',2)
            ylabel('Angular speed (deg/s)');
            title('Mean angular speed');
            ylim([0 max(mean_mean_ang_speed)+10]);
            xlim([0 4]);
            xticks(1:3);
            xticklabels({'Darkness','Low contrast','High contrast'});
            
            %save figure
            saveas(gcf,[path,'plots\closedLoopAngSpeedPerBlock.png']);
            
            %% Relationship between bump magnitude and velocity
            
            allBumpMag = [];
            for block = 1:length(blockLimits)
                allBumpMag = [allBumpMag,bump_mag{block}];
            end
            
            figure('Position',[200 200 1600 800]),
            subplot(1,4,1)
            scatter(data.vel_for_deg_ds(1:length(allBumpMag)),allBumpMag)
            xlim([0 max(data.vel_for_deg_ds)]);
            xlabel('Forward velocity (deg/s)');
            ylabel('Bump magnitude (max-min)');
            
            subplot(1,4,2)
            scatter(abs(data.vel_side_deg_ds(1:length(allBumpMag))),allBumpMag,[],[0.8 0.2 0.6])
            xlim([0 max(abs(data.vel_side_deg_ds))]);
            xlabel('Side speed (deg/s)');
            
            subplot(1,4,3)
            scatter(abs(data.vel_yaw_ds(1:length(allBumpMag))),allBumpMag,[],[0.4 0.2 0.8])
            xlim([0 max(abs(data.vel_yaw_ds))]);
            xlabel('Angular speed (deg/s)');
            
            subplot(1,4,4)
            scatter(data.total_mvt_ds(1:length(allBumpMag)),allBumpMag,[],'k')
            xlim([0 max(data.total_mvt_ds)]);
            xlabel('Total movement (deg/s)');
            
            %save figure
            saveas(gcf,[path,'plots\closedLoopBMvsVel.png']);
            
            %% Binning it
            
            figure('Position',[200 200 1400 600]),
            %Plot bump magnitude in time
            subplot(2,4,[1 3])
            plot(data.time(1:length(allBumpMag)),allBumpMag, 'k')
            ylabel('Bump magnitude');
            ylim([0 max(allBumpMag)+0.5]);
            set(gca,'xticklabel',{[]});
            
            %Plot total movement
            subplot(2,4,[5 7])
            plot(data.time(1:length(allBumpMag)),data.total_mvt_ds(1:length(allBumpMag)),'k')
            xlabel('Time (sec)');
            ylabel('Total movement (deg/s)');
            
            %Plot relationship between both parameters
            subplot(2,4,[4,8]);
            
            %Create vector with the contrast level for each timepoint
            all_contrast_levels = [];
            for block = 1:length(blockLimits)
                contrast_level{block} = repelem(Intensities(block),blockLimits{block}(2)+1-blockLimits{block}(1));
                all_contrast_levels = [all_contrast_levels,contrast_level{block}];
            end
            
            %Bin the bump magnitude data using both the total movement and the stimulus
            %type
            nbins = 20;
            maxBin = max(data.total_mvt_ds);
            binWidth = maxBin/nbins;
            mvtBins = [0:binWidth:maxBin];
            
            %create axes for plot
            mvtAxes = mvtBins - binWidth;
            mvtAxes = mvtAxes(2:end);
            mvtAxes(end) = mvtAxes(end-1)+binWidth;
            
            %Getting binned means
            for bin = 1:length(mvtBins)-1
                for contrast = 1:3
                    doubleBin(bin,contrast) = mean(allBumpMag((data.total_mvt_ds(1:length(allBumpMag)) > mvtBins(bin)) & (data.total_mvt_ds(1:length(allBumpMag)) < mvtBins(bin+1)) & (all_contrast_levels == contrast)));
                end
            end
            plot(mvtAxes,doubleBin,'-o')
            ylabel('Mean bump magnitude'); xlabel('Total movement (deg/s)');
            ylim([0 (max(max(doubleBin))+0.5)]);
            legend('Darkness','Low contrast','High contrast');
            
            %save figure
            saveas(gcf,[path,'plots\closedLoopBMvsVelPerContrastBinned.png']);
            
            %% Model bump magnitude as a function of contrast level and total movement
            
            %fill missing data in the total movement vector
            data.total_mvt_ds = fillmissing(data.total_mvt_ds,'linear');
            %zscore the total movement data
            zscored_total_mvt = zscore(data.total_mvt_ds);
            
            %Create table with the model's variables
            modelTable = table(all_contrast_levels',data.total_mvt_ds(1:length(allBumpMag))',zscored_total_mvt(1:length(allBumpMag))',data.time(1:length(allBumpMag)),allBumpMag','VariableNames',{'ContrastLevel','TotalMovement','ZscoredTotalMovement','Time','BumpMagnitude'});
            
            %Fit linear model using contrast level as a categorical variable
            mdl_BM = fitlm(modelTable,'BumpMagnitude~ContrastLevel+TotalMovement','CategoricalVars',1)
            %Using the zscored movement: this model is equal to the previous one, but
            %zscored data will be useful when pooling data across flies
            mdl_BMz = fitlm(modelTable,'BumpMagnitude~ContrastLevel+ZscoredTotalMovement','CategoricalVars',1)
            
            %% Relationship between bump width at half max and velocity
            
            allBumpWidth = [];
            for block = 1:length(blockLimits)
                allBumpWidth = [allBumpWidth,width_half_max{block}];
            end
            
            figure('Position',[200 200 1600 800]),
            subplot(1,4,1)
            scatter(data.vel_for_deg_ds(1:length(allBumpWidth)),allBumpWidth)
            xlim([0 max(data.vel_for_deg_ds)]);
            xlabel('Forward velocity (deg/s)');
            ylabel('Bump width at half max');
            
            subplot(1,4,2)
            scatter(abs(data.vel_side_deg_ds(1:length(allBumpWidth))),allBumpWidth,[],[0.8 0.2 0.6])
            xlim([0 max(abs(data.vel_side_deg_ds))]);
            xlabel('Side speed (deg/s)');
            
            subplot(1,4,3)
            scatter(abs(data.vel_yaw_ds(1:length(allBumpWidth))),allBumpWidth,[],[0.4 0.2 0.8])
            xlim([0 max(abs(data.vel_yaw_ds))]);
            xlabel('Angular speed (deg/s)');
            
            subplot(1,4,4)
            scatter(data.total_mvt_ds(1:length(allBumpWidth)),allBumpWidth,[],'k')
            xlim([0 max(data.total_mvt_ds)]);
            xlabel('Total movement (deg/s)');
            
            %save figure
            saveas(gcf,[path,'plots\closedLoopHWvsVel.png']);
            
            
            %% Binning it
            
            %Plot
            figure('Position',[200 200 1400 600]),
            %Plot bump magnitude in time
            subplot(2,4,[1 3])
            plot(data.time(1:length(allBumpWidth)),allBumpWidth, 'k')
            ylabel('Bump width at half max');
            ylim([0 8]);
            set(gca,'xticklabel',{[]});
            
            %Plot total movement
            subplot(2,4,[5 7])
            plot(data.time(1:length(allBumpWidth)),data.total_mvt_ds(1:length(allBumpWidth)),'k')
            xlabel('Time (sec)');
            ylabel('Total movement (deg/s)');
            
            %Plot relationship between both parameters
            subplot(2,4,[4,8]);
            
            %getting binned means
            for bin = 1:length(mvtBins)-1
                for contrast = 1:3
                    doubleBin(bin,contrast) = mean(allBumpWidth((data.total_mvt_ds(1:length(allBumpWidth)) > mvtBins(bin)) & (data.total_mvt_ds(1:length(allBumpWidth)) < mvtBins(bin+1)) & (all_contrast_levels == contrast)));
                end
            end
            plot(mvtAxes,doubleBin,'-o')
            ylabel('Mean bump width at half max'); xlabel('Total movement (deg/s)');
            ylim([0 (max(max(doubleBin))+0.5)]);
            legend('Darkness','Low contrast','High contrast');
            
            %save figure
            saveas(gcf,[path,'plots\closedLoopHWvsVelPerContrastBinned.png']);
            
            %% Model bump width at half max as a function of contrastLevel and total movement
            
            %create table with the model's variables
            modelTable_HW = table(all_contrast_levels',data.total_mvt_ds(1:length(allBumpWidth))',zscored_total_mvt(1:length(allBumpWidth))',data.time(1:length(allBumpWidth)),allBumpWidth','VariableNames',{'ContrastLevel','TotalMovement','ZscoredTotalMovement','Time','BumpWidth'});
            
            %fit linear model using contrast level as a categorical variable
            mdl_HW = fitlm(modelTable_HW,'BumpWidth~ContrastLevel+TotalMovement','CategoricalVars',1)
            mdl_HW2 = fitlm(modelTable_HW,'BumpWidth~ContrastLevel+ZscoredTotalMovement','CategoricalVars',1)
            
            %% Heading variation per stimulus
            
            figure('Position',[200 200 1000 800]),
            subplot(2,1,1)
            for block = 1:length(blockLimits)
                [~, heading_var(block)] = circ_std(deg2rad(heading(blockLimits{block}(1):blockLimits{block}(2))),[],[],1);
                plot(block,heading_var(block),'ko','MarkerFaceColor',color_gradient{Intensities(block)},'MarkerSize',8)
                hold on
            end
            xlim([0 num_subplots+1]);
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
            saveas(gcf,[path,'plots\closedLoopHeadingVariation.png']);
            
            
            %% Bump stability analysis
            
            %we're going to use bump stability as a proxy for intended menotaxis
            %Let's first pool all the data with a certain contrast level
            bump_pos = data.dff_pva;
            for contrast = 1:3
                per_contrast_data{contrast,1} = bump_pos(all_contrast_levels==contrast);
                per_contrast_data{contrast,2} = allBumpMag(all_contrast_levels==contrast);
                per_contrast_data{contrast,3} = allBumpWidth(all_contrast_levels==contrast);
            end
            
            %define the possible bins
            bins = [5,10,15,20,30,40,50,100];
            
            %compute relationship between bump stability and bump parameters for each
            %contrast and temporal window
            for contrast = 1:3
                for bin_amount = 1:length(bins)
                    
                    %reset variables
                    clear binned_bump_mag
                    clear binned_half_width
                    clear corr_coefficient_bm
                    clear corr_coefficient_hw
                    clear pvalue_bm
                    clear pvalue_hw
                    clear bump_var
                    
                    %change bin number
                    bin_number = bins(bin_amount);
                    bin_width = floor(length(per_contrast_data{contrast,1})/bin_number);
                    
                    %get bin limits
                    bin_limits{1} = [1,bin_width+1];
                    for bin = 2:bin_number-1
                        bin_limits{bin} = [2+bin_width*(bin-1),1+bin_width*(bin-1)+bin_width];
                    end
                    bin_limits{bin_number} = [length(per_contrast_data{contrast,1})-bin_width,length(per_contrast_data{contrast,1})];
                    
                    %get circular standard deviation of bump position
                    for bin = 1:bin_number
                        [~, bump_var(bin)] = circ_std(per_contrast_data{contrast,1}(bin_limits{bin}(1):bin_limits{bin}(2)));
                    end
                    
                    %get bump magnitude per bin
                    for bin = 1:bin_number
                        binned_bump_mag(bin) = median(per_contrast_data{contrast,2}(bin_limits{bin}(1):bin_limits{bin}(2)));
                    end
                    
                    %get bump width at half max per bin
                    for bin = 1:bin_number
                        binned_half_width(bin) = median(per_contrast_data{contrast,3}(bin_limits{bin}(1):bin_limits{bin}(2)));
                    end
                    
                    %get correlation between variables
                    [corr_coefficient_bm,pvalue_bm] = corrcoef(bump_var,binned_bump_mag);
                    correlations_bm{contrast}(bin_amount) = corr_coefficient_bm(1,2);
                    pval_bm{contrast}(bin_amount) = pvalue_bm(1,2);
                    
                    [corr_coefficient_hw,pvalue_hw] = corrcoef(bump_var,binned_half_width);
                    correlations_hw{contrast}(bin_amount) = corr_coefficient_hw(1,2);
                    pval_hw{contrast}(bin_amount) = pvalue_hw(1,2);
                    
                end
            end
            
            %Plot the relationships
            figure('Position',[200 200 1200 800]),
            for contrast = 1:3
                subplot(3,2,contrast*2-1)
                plot(bins(pval_bm{contrast}<0.05),correlations_bm{contrast}(pval_bm{contrast}<0.05),'ro')
                hold on
                plot(bins(pval_bm{contrast}>=0.05),correlations_bm{contrast}(pval_bm{contrast}>=0.05),'ko')
                line([0 100],[0 0],'color','b');
                ylim([-1 1]);
                xlabel('Bin number');
                ylabel('Correlation coefficient');
                title('Relationship between bump variation and bump magnitude');
                
                subplot(3,2,contrast*2)
                plot(bins(pval_hw{contrast}<0.05),correlations_hw{contrast}(pval_hw{contrast}<0.05),'ro')
                hold on
                plot(bins(pval_hw{contrast}>=0.05),correlations_hw{contrast}(pval_hw{contrast}>=0.05),'ko')
                line([0 100],[0 0],'color','b');
                ylim([-1 1]);
                xlabel('Bin number');
                ylabel('Correlation coefficient');
                title('Relationship between bump variation and bump width at half max');
            end
            
            
            %save figure
            saveas(gcf,[path,'plots\bumpStabilityAnalysis.png']);
            
            %% Repeat bump stability analysis using a threshold for movement
            
            %find threshold using total movement distribution
            
            figure
            h = histogram(data.total_mvt_ds,30,'Normalization','probability');
            kernelfit = fitdist(data.total_mvt_ds','Kernel','BandWidth',8);
            x = 0:max(data.total_mvt_ds)/30:max(data.total_mvt_ds-1);
            data_pdf = pdf(kernelfit,x);
            hold on
            plot(x,data_pdf*10,'k-','LineWidth',2)
            
            %make sure the histogram starts at 0
            h.BinLimits = [0 h.BinLimits(2)]
            h.NumBins = 30
            %get total movement value
            if contains(path, '20201201_60D05_7f_fly2')
                density_min = 4;
            elseif contains(path, '20201208_60D05_7f_fly2')
                density_min = 2;
            else
                %find local minima
                density_min = find(islocalmin(data_pdf));
            end
            mvt_thresh = h.BinEdges(density_min(1));
            
            %add to plot to check
            plot(mvt_thresh,data_pdf(density_min(1))*10,'ro','markerfacecolor','r')
            
            %cut the data to match the length up to the last block limit
            total_mvt = data.total_mvt_ds(1:length(all_contrast_levels));
            for contrast = 1:3
                per_contrast_data_thresh{contrast,1} = bump_pos(all_contrast_levels==contrast & total_mvt > mvt_thresh);
                per_contrast_data_thresh{contrast,2} = allBumpMag(all_contrast_levels==contrast & total_mvt > mvt_thresh);
                per_contrast_data_thresh{contrast,3} = allBumpWidth(all_contrast_levels==contrast & total_mvt > mvt_thresh);
            end
            
            %compute relationship between bump stability and bump parameters for each
            %contrast and temporal window
            for contrast = 1:3
                for bin_amount = 1:length(bins)
                    
                    %reset variables
                    clear corr_coefficient_bm
                    clear corr_coefficient_hw
                    clear pvalue_bm
                    clear pvalue_hw
                    
                    %change bin number
                    bin_number = bins(bin_amount);
                    bin_width = floor(length(per_contrast_data_thresh{contrast,1})/bin_number);
                    
                    %get bin limits
                    bin_limits{1} = [1,bin_width+1];
                    for bin = 2:bin_number-1
                        bin_limits{bin} = [2+bin_width*(bin-1),1+bin_width*(bin-1)+bin_width];
                    end
                    bin_limits{bin_number} = [length(per_contrast_data_thresh{contrast,1})-bin_width,length(per_contrast_data_thresh{contrast,1})];
                    
                    %get circular standard deviation of bump position
                    for bin = 1:bin_number
                        [~, bump_var_per_contrast{contrast,bin_amount}(bin)] = circ_std(per_contrast_data_thresh{contrast,1}(bin_limits{bin}(1):bin_limits{bin}(2)));
                    end
                    
                    %get bump magnitude per bin
                    for bin = 1:bin_number
                        binned_bump_mag_per_contrast{contrast,bin_amount}(bin) = median(per_contrast_data_thresh{contrast,2}(bin_limits{bin}(1):bin_limits{bin}(2)));
                    end
                    
                    %get bump width at half max per bin
                    for bin = 1:bin_number
                        binned_half_width_per_contrast{contrast,bin_amount}(bin) = median(per_contrast_data_thresh{contrast,3}(bin_limits{bin}(1):bin_limits{bin}(2)));
                    end
                    
                    %get correlation between variables
                    [corr_coefficient_bm,pvalue_bm] = corrcoef(bump_var_per_contrast{contrast,bin_amount},binned_bump_mag_per_contrast{contrast,bin_amount});
                    correlations_bm_thresh{contrast}(bin_amount) = corr_coefficient_bm(1,2);
                    pval_bm_thresh{contrast}(bin_amount) = pvalue_bm(1,2);
                    
                    [corr_coefficient_hw,pvalue_hw] = corrcoef(bump_var_per_contrast{contrast,bin_amount},binned_half_width_per_contrast{contrast,bin_amount});
                    correlations_hw_thresh{contrast}(bin_amount) = corr_coefficient_hw(1,2);
                    pval_hw_thresh{contrast}(bin_amount) = pvalue_hw(1,2);
                    
                end
            end
            
            %Plot the relationships
            figure('Position',[200 200 1200 800]),
            for contrast = 1:3
                subplot(3,2,contrast*2-1)
                plot(bins(pval_bm_thresh{contrast}<0.05),correlations_bm_thresh{contrast}(pval_bm_thresh{contrast}<0.05),'ro')
                hold on
                plot(bins(pval_bm_thresh{contrast}>=0.05),correlations_bm_thresh{contrast}(pval_bm_thresh{contrast}>=0.05),'ko')
                line([0 100],[0 0],'color','b');
                ylim([-1 1]);
                xlabel('Bin number');
                ylabel('Correlation coefficient');
                title('Relationship between bump variation and bump magnitude');
                
                subplot(3,2,contrast*2)
                plot(bins(pval_hw_thresh{contrast}<0.05),correlations_hw_thresh{contrast}(pval_hw_thresh{contrast}<0.05),'ro')
                hold on
                plot(bins(pval_hw_thresh{contrast}>=0.05),correlations_hw_thresh{contrast}(pval_hw_thresh{contrast}>=0.05),'ko')
                line([0 100],[0 0],'color','b');
                ylim([-1 1]);
                xlabel('Bin number');
                ylabel('Correlation coefficient');
                title('Relationship between bump variation and bump width at half max');
            end
            
            %save figure
            saveas(gcf,[path,'plots\bumpStabilityAnalysisAboveThresh.png']);
            
            %% Repeat analysis binning data from all blocks
            
            pooled_data_thresh{1} = bump_pos(total_mvt > mvt_thresh);
            pooled_data_thresh{2} = allBumpMag(total_mvt > mvt_thresh);
            pooled_data_thresh{3} = allBumpWidth(total_mvt > mvt_thresh);
            
            clear bump_var
            clear binned_bump_mag
            clear binned_half_width
            
            %compute relationship between bump stability and bump parameters for each
            %contrast and temporal window
            for bin_amount = 1:length(bins)
                
                %reset variables
                clear corr_coefficient_bm
                clear corr_coefficient_hw
                clear pvalue_bm
                clear pvalue_hw
                
                %change bin number
                bin_number = bins(bin_amount);
                bin_width = floor(length(pooled_data_thresh{1})/bin_number);
                
                %get bin limits
                bin_limits{1} = [1,bin_width+1];
                for bin = 2:bin_number-1
                    bin_limits{bin} = [2+bin_width*(bin-1),1+bin_width*(bin-1)+bin_width];
                end
                bin_limits{bin_number} = [length(pooled_data_thresh{1})-bin_width,length(pooled_data_thresh{1})];
                
                %get circular standard deviation of bump position
                for bin = 1:bin_number
                    [~, bump_var{bin_amount}(bin)] = circ_std(pooled_data_thresh{1}(bin_limits{bin}(1):bin_limits{bin}(2)));
                end
                
                %get bump magnitude per bin
                for bin = 1:bin_number
                    binned_bump_mag{bin_amount}(bin) = median(pooled_data_thresh{2}(bin_limits{bin}(1):bin_limits{bin}(2)));
                end
                
                %get bump width at half max per bin
                for bin = 1:bin_number
                    binned_half_width{bin_amount}(bin) = median(pooled_data_thresh{3}(bin_limits{bin}(1):bin_limits{bin}(2)));
                end
                
                %get correlation between variables
                [corr_coefficient_bm,pvalue_bm] = corrcoef(bump_var{bin_amount},binned_bump_mag{bin_amount});
                correlations_bm_thresh_pooled(bin_amount) = corr_coefficient_bm(1,2);
                pval_bm_thresh_pooled(bin_amount) = pvalue_bm(1,2);
                
                [corr_coefficient_hw,pvalue_hw] = corrcoef(bump_var{bin_amount},binned_half_width{bin_amount});
                correlations_hw_thresh_pooled(bin_amount) = corr_coefficient_hw(1,2);
                pval_hw_thresh_pooled(bin_amount) = pvalue_hw(1,2);
                
            end
            
            %Plot the relationships
            figure('Position',[200 200 1200 800]),
            subplot(1,2,1)
            plot(bins(pval_bm_thresh_pooled<0.05),correlations_bm_thresh_pooled(pval_bm_thresh_pooled<0.05),'ro')
            hold on
            plot(bins(pval_bm_thresh_pooled>=0.05),correlations_bm_thresh_pooled(pval_bm_thresh_pooled>=0.05),'ko')
            line([0 100],[0 0],'color','b');
            ylim([-1 1]);
            xlabel('Bin number');
            ylabel('Correlation coefficient');
            title('Relationship between bump variation and bump magnitude');
            
            subplot(1,2,2)
            plot(bins(pval_hw_thresh_pooled<0.05),correlations_hw_thresh_pooled(pval_hw_thresh_pooled<0.05),'ro')
            hold on
            plot(bins(pval_hw_thresh_pooled>=0.05),correlations_hw_thresh_pooled(pval_hw_thresh_pooled>=0.05),'ko')
            line([0 100],[0 0],'color','b');
            ylim([-1 1]);
            xlabel('Bin number');
            ylabel('Correlation coefficient');
            title('Relationship between bump variation and bump width at half max');
            
            %save figure
            saveas(gcf,[path,'plots\pooledBumpStabilityAnalysisAboveThresh.png']);
            
            %% Run bump stability model
            
            for bin_amount = 1:length(bins)
                bs_table{bin_amount} = table(bump_var{1,bin_amount}',binned_bump_mag{1,bin_amount}',binned_half_width{1,bin_amount}','VariableNames',{'bump_var','bump_mag','bump_width'});
                %fitlm(bs_table{bin_amount},'bump_var~bump_mag+bump_width-1')
            end
            
            
            %% Save useful data
            
            save([path,'\summary_data.mat'],'summary_data','mean_reference_offset','modelTable','modelTable_HW','correlations_hw','correlations_bm','pval_hw','pval_bm','correlations_hw_thresh','correlations_bm_thresh','pval_hw_thresh','pval_bm_thresh','correlations_hw_thresh_pooled','correlations_bm_thresh_pooled','pval_hw_thresh_pooled','pval_bm_thresh_pooled','bs_table','bump_var_per_contrast','binned_bump_mag_per_contrast','binned_half_width_per_contrast')
            
            close all; clc;
            
        else
        end
    end
    
end
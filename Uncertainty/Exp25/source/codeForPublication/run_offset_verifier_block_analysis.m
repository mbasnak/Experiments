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
            
        end
    end
    
    close all;
end
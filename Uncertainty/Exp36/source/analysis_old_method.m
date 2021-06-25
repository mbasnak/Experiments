%Code to analyze experiment 36, in which the flies receive wind and visual
%cues that might be high or low reliability, and might have the same or
%different gain.

%clean up space
clear all; close all;

%get data directory
[path] = uigetdir('Z:\Wilson Lab\Mel\Experiments\Uncertainty\Exp36\data');

%% For each file in directory, import and plot, checking gain and cue
%reliability

fileNames = dir([path,'\analysis\']);

for file = 1:length(fileNames)
    if (contains(fileNames(file).name,'analysis') & ~contains(fileNames(file).name,'continuous'))
        %load data
        load(fullfile(fileNames(file).folder,fileNames(file).name))
        
        %load run_obj
        run_obj_dir = [path,'\ball\runobj\'];
        expression = ['*_sid_',num2str(data.sid),'_runobj.mat'];
        run_obj_file_name = dir(fullfile(run_obj_dir,expression));
        load(fullfile(run_obj_file_name.folder,run_obj_file_name.name))
        
        
        figure('Position',[100 100 1200 800]),
        subplot(3,6,[2 5])
        dff = data.dff_matrix;
        imagesc(dff)
        colormap(flipud(gray))
        title('EPG activity');
        
        subplot(3,6,[8 11])
        bump_pos = wrapTo180(rad2deg(data.phase));
        %Remove wrapped lines to plot
        [x_out_bump,bump_to_plot] = removeWrappedLines(data.time,bump_pos');
        plot(x_out_bump,bump_to_plot,'LineWidth',1.5)
        hold on
        panels = wrapTo180(data.panel_angle);
        [x_out_panels,panels_to_plot] = removeWrappedLines(data.time,panels);
        %Plot according to the bar's reliability
        if run_obj.pattern_number == 57
            plot(x_out_panels,panels_to_plot,'LineWidth',1.5,'color',[ 0.5 0.8 0.9])
        elseif run_obj.pattern_number == 56
            plot(x_out_panels,panels_to_plot,'LineWidth',1.5,'color',[ 0 0 0.6])           
        end
        wind = wrapTo180(rad2deg(data.motor_pos));
        [x_out_wind,wind_to_plot] = removeWrappedLines(data.time,wind');
        %Plot according to the wind's reliability
        if run_obj.airflow.Value == 0.15
            plot(x_out_wind,wind_to_plot,'LineWidth',1.5,'color',[.5 .5 .5])
        elseif run_obj.airflow.Value == 0.05
            plot(x_out_wind,wind_to_plot,'LineWidth',1.5,'color',[0 0 0])            
        end
        title('Stimuli and fly position');
        legend('EPG phase','Bar position','Wind position');
        ylim([-180 180]);
        
        subplot(3,6,[14 17])
        panels_offset = wrapTo180(rad2deg(circ_dist(data.phase',deg2rad(data.panel_angle))));
        [x_out_panels_offset,panels_offset_to_plot] = removeWrappedLines(data.time,panels_offset);
        plot(x_out_panels_offset,panels_offset_to_plot,'LineWidth',1.5)
        wind_offset = wrapTo180(rad2deg(circ_dist(data.phase,data.motor_pos)));
        [x_out_wind_offset,wind_offset_to_plot] = removeWrappedLines(data.time,wind_offset');
        hold on
        plot(x_out_wind_offset,wind_offset_to_plot,'LineWidth',1.5)
        title('Panels and wind offset');
        legend('Panels offset','Wind offset');
        ylim([-180 180]);
        
        subplot(3,6,13)
        polarhistogram(deg2rad(panels_offset))
        title('Panels offset');
        
        subplot(3,6,18)
        polarhistogram(deg2rad(wind_offset))
        title('Wind offset');
        
        
        %Define the cue reliability
        if run_obj.pattern_number == 57
            bar_reliability = 'reliable';
        else
            bar_reliability = 'unreliable';
        end
        if run_obj.airflow.Value == 0.15
            wind_reliability = 'reliable';
        else
            wind_reliability = 'unreliable';
        end
        
        
        suptitle([bar_reliability,' visual stimulus, ',wind_reliability,' wind stimulus']);
        
        saveas(gcf,[path,'\analysis\plots\heatmap_and_offset_',bar_reliability,'_bar_gain',num2str(run_obj.gain_panels),wind_reliability,'_wind_gain',num2str(run_obj.gain_wind),'.png']);
        
        
        %% Correlation between cue positions
        
        figure,
        scatter(panels,wind)
        xlabel('Bar position');
        ylabel('Wind position');
        
        saveas(gcf,[path,'\analysis\plots\cue_pos_correlation_',bar_reliability,'_bar_gain',num2str(run_obj.gain_panels),wind_reliability,'_wind_gain',num2str(run_obj.gain_wind),'.png']);
        
    end
end


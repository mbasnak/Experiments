%code to analyze experiment 40 (offset control)

close all; clear all;

%get directory
[path] = uigetdir('Z:\Wilson Lab\Mel\Experiments\Uncertainty\Exp40\data');
fileNames = dir([path,'\analysis']);

%get session information
session_info = load([path,'\analysis\session_info.mat']);

%% Make dir to save plots

mkdir([path,'\analysis'],'\plots')

%% 

%load files
for file = 1:length(fileNames)
    
    if contains(fileNames(file).name,'continuous')
        
        %load the data
        load(fullfile(fileNames(file).folder,fileNames(file).name))
               
        %Plot full trial
        figure('Position',[100 100 1800 1000]),
        subplot(5,1,1)
        dff = continuous_data.dff_matrix';
        imagesc(flip(dff))
        colormap(flipud(gray))
        title('EPG activity');
        set(gca,'xticklabel',{[]});
        set(gca,'yticklabel',{[]});
        
        subplot(5,1,2)
        bump_pos = wrapTo180(rad2deg(continuous_data.bump_pos));
        %Remove wrapped lines to plot
        [x_out_bump,bump_to_plot] = removeWrappedLines(continuous_data.time,bump_pos');
        plot(x_out_bump,bump_to_plot,'LineWidth',1.5)
        hold on
        %add fly position
        fly_pos = wrapTo180(-continuous_data.heading_deg);
        [x_out_fly,fly_to_plot] = removeWrappedLines(continuous_data.time,fly_pos);
        plot(x_out_fly,fly_to_plot,'.','LineWidth',1.5)
        title('Bump and fly position');
        legend('Bump estimate','Fly position')
        ylim([-180 180]);
        xlim([0 x_out_bump(end)]);
        set(gca,'xticklabel',{[]})
        
        subplot(5,1,3)
        offset = wrapTo180(rad2deg(circ_dist(continuous_data.bump_pos',-continuous_data.heading)));
        [x_out_offset,offset_to_plot] = removeWrappedLines(continuous_data.time,offset);
        plot(x_out_offset,offset_to_plot,'LineWidth',1.5,'color',[0.8500 0.3250 0.0980])
        title('Offset');
        ylim([-180 180]);
        xlim([0 x_out_offset(end-1)]);
        set(gca,'xticklabel',{[]})
        
        subplot(5,1,4)
        plot(continuous_data.time(continuous_data.adj_rs>=0.5),continuous_data.bump_magnitude(continuous_data.adj_rs>=0.5),'r.','handleVisibility','off')
        hold on
        plot(continuous_data.time(continuous_data.adj_rs<0.5),continuous_data.bump_magnitude(continuous_data.adj_rs<0.5),'k.','handleVisibility','off')
        title('Bump magnitude')
        xlim([0 continuous_data.time(end)]);
        set(gca,'xticklabel',{[]})
        
        subplot(5,1,5)
        plot(continuous_data.time(continuous_data.adj_rs>=0.5),continuous_data.bump_width(continuous_data.adj_rs>=0.5),'r.','handleVisibility','off')
        hold on
        plot(continuous_data.time(continuous_data.adj_rs<0.5),continuous_data.bump_width(continuous_data.adj_rs<0.5),'k.','handleVisibility','off')
        xlabel('Time (sec)');
        title('Bump width');
        xlim([0 continuous_data.time(end)]);
        
        %store offset
        %determine type of trial
        sid = str2num(fileNames(file).name(25));
        
        if sid == session_info.session_info.empty
            trial = 'empty';
            offset_empty = offset;
        elseif sid == session_info.session_info.bar
            trial = 'bar';
            offset_bar = offset;
        else
            trial = 'wind';
            offset_wind = offset;
        end
        

    end
       
end

%% Plot offset distribution

%Plot offset distribution
figure,
subplot(1,3,1)
polarhistogram(deg2rad(offset_empty))
title('Empty trial');

subplot(1,3,2)
polarhistogram(deg2rad(offset_wind))
title('Wind trial');

subplot(1,3,3)
polarhistogram(deg2rad(offset_bar))
title('Bar trial');


%% Get and plot offset variability 

[~, offset_var_empty] = circ_std(deg2rad(offset_empty));
[~, offset_var_wind] = circ_std(deg2rad(offset_wind));
[~, offset_var_bar] = circ_std(deg2rad(offset_bar));

all_offset_var = [offset_var_empty,offset_var_wind,offset_var_bar];

figure,
plot(all_offset_var,'-ko')
xlim([0 4]);
ylim([0 3]);
xticks([1:3]);
xticklabels({'Empty trial','Wind trial','Bar trial'});
ylabel('Offset variability (rad)');

saveas(gcf,[path,'\analysis\plots\offset_var.png'])

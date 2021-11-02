%analysis of the initial offset stabilizing block


%% Load data

clear all; close all;

%Get the pre-processed data
[path] = uigetdir('Z:\Wilson Lab\Mel\Experiments\Uncertainty\Exp37\data\');

% Import sessions information
load([path,'\analysis\sessions_info.mat'])
load([path,'\analysis\continuous_analysis_sid_',num2str(sessions.offset_stabilizer),'_tid_0.mat'])


%% Make directory to save plots

%Move to the analysis folder
cd([path,'\analysis\'])
%List the contents
contents = dir();
%if there isn't a 'plots' folder already, create one
if (contains([contents.name],'plots') == 0)
   mkdir(path,'\analysis\plots'); 
end
%List the contents of the 'plots' folder
cd([path,'\analysis\plots\'])

%% Set colormap

folderNames = dir(path(1:65));
flyNames = struct();
for folder = 1:length(folderNames)
    if (contains(folderNames(folder).name,'60D05') & ~contains(folderNames(folder).name,'txt'))
        flyNames(folder).name = folderNames(folder).name;
    end
end
%Remove empty rows
flyNames = flyNames(~cellfun(@isempty,{flyNames.name}));

%Assign fly number
for fly = 1:length(flyNames)
    if strcmp(flyNames(fly).name,path(66:end))
    %if strcmp(flyNames(fly).name,path(71:end))
        fly_ID = fly;
    end
end

colors_for_plots = [0.2 0.8 0.8 ; 1 0.5 0; 0 0.5 1;...
    0 0.6 0.3;  1 0.2 0.2; 0.9290 0.6940 0.1250;...
    0.6350 0.0780 0.1840; 0.4660 0.6740 0.1880;...
    0 0.4470 0.7410; 0.75, 0.1, 0.75; 0.75, 0.75, 0];

fly_color = colors_for_plots(fly_ID,:);


%% Plot activity heatmap with bump parameters

figure('Position',[100 100 1200 800]),
subplot(5,1,1)
dff = continuous_data.dff_matrix';
imagesc(flip(dff))
colormap(flipud(gray))
title('EPG activity');
set(gca,'xticklabel',{[]})
set(gca,'yticklabel',{[]})

subplot(5,1,2)
bump_pos = wrapTo180(rad2deg(continuous_data.bump_pos));
%Remove wrapped lines to plot
[x_out_bump,bump_to_plot] = removeWrappedLines(continuous_data.time,bump_pos');
plot(x_out_bump,bump_to_plot,'LineWidth',1.5)
hold on
heading = wrapTo180(-continuous_data.heading_deg);
[x_out_heading,heading_to_plot] = removeWrappedLines(continuous_data.time,heading);
plot(x_out_heading,heading_to_plot,'LineWidth',1.5)
title('Bump and fly position');
legend('Bump estimate','Fly position','Location','best')
ylim([-180 180]);
set(gca,'xticklabel',{[]})

subplot(5,1,3)
offset = wrapTo180(rad2deg(circ_dist(continuous_data.bump_pos',-continuous_data.heading)));
[x_out_offset,offset_to_plot] = removeWrappedLines(continuous_data.time,offset);
plot(x_out_offset,offset_to_plot,'LineWidth',1.5,'color','k')
title('Offset')
ylim([-180 180]);
set(gca,'xticklabel',{[]})

subplot(5,1,4)
plot(continuous_data.time(continuous_data.adj_rs>=0.5),continuous_data.bump_magnitude(continuous_data.adj_rs>=0.5),'r.')
hold on
plot(continuous_data.time(continuous_data.adj_rs<0.5),continuous_data.bump_magnitude(continuous_data.adj_rs<0.5),'k.')
title('Bump magnitude')

subplot(5,1,5)
plot(continuous_data.time(continuous_data.adj_rs>=0.5),continuous_data.bump_width(continuous_data.adj_rs>=0.5),'r.')
hold on
plot(continuous_data.time(continuous_data.adj_rs<0.5),continuous_data.bump_width(continuous_data.adj_rs<0.5),'k.')
title('Bump width')

saveas(gcf,[path,'\analysis\plots\offset_stabilizing_block.png']);

%%
clear all; close all;
%Code to analyze the data for the bar contrast block across flies

%% Load the data
clear all; close all;
%Get directory
dname = uigetdir('Z:\Wilson Lab\Mel\Experiments\Uncertainty\Exp25\data\Experimental');
%Find data folders
dir_contents = dir(dname);
for content = 1:length(dir_contents)
    if contains(dir_contents(content).name,'60D05')
        %Check for 'bar_contrast' data files.
         expression = ['bar_contrast_data' '*.mat'];
         data_dir = [dir_contents(content).folder, '\', dir_contents(content).name, '\analysis\'];
         analysis_file = dir(fullfile(data_dir, expression));
         %Load as field in struct
         if length(analysis_file) == 1
            data(content) = load(fullfile(analysis_file.folder, analysis_file.name));
         else
            data(content) = load(fullfile(analysis_file(1).folder, analysis_file(1).name));
            extra(content) = load(fullfile(analysis_file(2).folder, analysis_file(2).name));
         end
    end
end

%Remove empty rows
empty_elems = arrayfun(@(s) all(structfun(@isempty,s)), data);
data = data(~empty_elems);

if exist('extra','var')
    empty_elems2 = arrayfun(@(s) all(structfun(@isempty,s)), extra);
    extra = extra(~empty_elems2);
    data = [data,extra];
end

%% Analyze change in offset

figure('Position',[200, 200, 1200, 800]),
subplot(1,2,1)
for fly = 1:length(data)
    plot(data(fly).orderedOffsetChange,'o')
    hold on
end
legend({'Fly1','Fly2','Fly3'});
ylim([-180 180]);
yline(0,'--r','HandleVisibility','off');
xlim([0 7]);
xticks([1:6]);
xlabel('Contrast level');
ylabel('Offset change around bar jump (deg)');

%Create table to fit regression
pre_jump_offset = [];
post_jump_offset = [];
for fly = 1: length(data)
    pre_jump_offset = [pre_jump_offset; data(fly).ordered_pre_jump_offset(2:end)];
    post_jump_offset = [post_jump_offset; data(fly).ordered_post_jump_offset(2:end)];
end
offsetChange = array2table([pre_jump_offset,post_jump_offset], 'VariableNames', {'pre_jump_offset', 'post_jump_offset'});
%Fit a model to all contrasts excluding darkness 
offsetChangeModel = fitlm(offsetChange);

subplot(1,2,2)
for fly = 1:length(data)
    plot(data(fly).ordered_pre_jump_offset(2:end),data(fly).ordered_post_jump_offset(2:end),'o')
    hold on
end
xlim([-180 180]); ylim([-180 180]);
xlabel('Offset before cue jump'); ylabel('Offset after cue jump');
%add regression line
hold on
mock_offsets = -180:10:180;
regression_line = offsetChangeModel.Coefficients.Estimate(1) + offsetChangeModel.Coefficients.Estimate(2)*mock_offsets;
plot(mock_offsets,regression_line,'r')
%Add R and p-value as text
text(100,120,['R =',num2str(round(offsetChangeModel.Rsquared.Ordinary,2,'significant'))]);
text(100,140,['p-val =',num2str(round(offsetChangeModel.Coefficients.pValue(2),2,'significant'))]);

%Save figure:
saveas(gcf,[dname,'\globalPlots\OffsetChange.png']);

%% Offset variation vs contrast

figure,
for fly = 1:length(data)
    plot(data(fly).offset_var,'color',[.6 .6 .6])
    hold on
end
offset_var = [];
for fly = 1: length(data)
    offset_var = [offset_var; data(fly).offset_var];
end
plot(mean(offset_var),'k','LineWidth',2)
xlabel('Contrast level');
ylabel('Offset variation (circular std)');
xlim([0.5 6.5]);
xticks(1:6);

%Save figure:
saveas(gcf,[dname,'\globalPlots\OffsetVariation.png']);

%% Heading variation vs contrast

figure,
for fly = 1:length(data)
    plot(data(fly).heading_var,'color',[.6 .6 .6])
    hold on
end
heading_var = [];
for fly = 1: length(data)
    heading_var = [heading_var; data(fly).heading_var];
end
plot(mean(heading_var),'k','LineWidth',2)
xlabel('Contrast level');
ylabel('Heading variation (circular std)');
xlim([0.5 6.5]);
xticks(1:6);

%Save figure:
saveas(gcf,[dname,'\globalPlots\HeadingVariation.png']);

%% Bump mag vs contrast

figure,
subplot(2,1,1)
for fly = 1:length(data)
    plot(data(fly).median_zscore_bump_mag,'color',[.6 .6 .6])
    hold on
end
zscore_bump_mag = [];
for fly = 1: length(data)
    zscore_bump_mag = [zscore_bump_mag; data(fly).median_zscore_bump_mag];
end
plot(mean(zscore_bump_mag),'k','LineWidth',2)
xlim([0.5 6.5]);
xticks(1:6);
ylabel('Median z-scored bump');

subplot(2,1,2)
for fly = 1:length(data)
    plot(data(fly).ordered_median_bump_mag,'color',[.6 .6 .6])
    hold on
end
vonMises_bump_mag = [];
for fly = 1: length(data)
    vonMises_bump_mag = [vonMises_bump_mag; data(fly).ordered_median_bump_mag];
end
plot(mean(vonMises_bump_mag),'k','LineWidth',2)
xlim([0.5 6.5]);
xticks(1:6);
ylabel('Median von Mises bump');
xlabel('Contrast level');

%Save figure:
saveas(gcf,[dname,'\globalPlots\BumpVsContrast.png']);

%% z-scored bump mag vs fwd vel and ang speed

%% changes in bump magnitude


figure,
for fly = 1:length(data)
    plot(data(fly).changes_in_bump_mag,'color',[.6 .6 .6])
    hold on
end
changes_in_bump_mag = [];
for fly = 1: length(data)
    changes_in_bump_mag = [changes_in_bump_mag; data(fly).changes_in_bump_mag];
end
plot(mean(changes_in_bump_mag),'k','LineWidth',2)
xlabel('Contrast level');
ylabel('Change in bump magnitude');
xlim([0.5 6.5]);
xticks(1:6);
yline(0,'-r');

%Save figure:
saveas(gcf,[dname,'\globalPlots\ChangeInBM.png']);
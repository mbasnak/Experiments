
clear all; close all

%% Load data

folderNames = dir('Z:\Wilson Lab\Mel\Experiments\Uncertainty\Exp30\data\changing_contrast');

offset = [];
longer_offset = [];
first_frame = [];

zscored_BM = [];
zscored_BW = [];
for folder = 1:length(folderNames)
   if contains(folderNames(folder).name,'60D05') == 1
       load(fullfile(folderNames(folder).folder,folderNames(folder).name,'analysis\offset_change.mat'))
       offset = [offset;offset_AJ];
       longer_offset = [longer_offset;longer_offset_AJ];
       first_frame = [first_frame;zscore(first_return_frame')];
       zscored_BM = [zscored_BM,zscore(mean_bm_pre_jump)];
       zscored_BW = [zscored_BW,zscore(mean_bw_pre_jump)];
   end
end

path = 'Z:\Wilson Lab\Mel\Experiments\Uncertainty\Exp30\data\changing_contrast\group_plots';


%% sort by bump paramters pre cue jump

%get indexes for sorted BM
[sorted_BM,I] = sort(zscored_BM);
[sorted_BW,Iw] = sort(zscored_BW);

%% Look at offset instead of change in offset, but shift for it to start at 0

shifted_offset = wrapTo180(offset-offset(:,1));
short_time_AJ = time_around_jump(376:426);

figure('Position',[100 100 1000 800]),
imagesc(shifted_offset)
colorbar
colormap(flipud(gray))
hold on
xline(26.5,'linewidth',2,'color','r');
xt = get(gca, 'XTick');                                            
xtlbl = linspace(min(short_time_AJ), max(short_time_AJ), numel(xt));                    
set(gca,'XTick',xt-2.5, 'XTickLabel',round(xtlbl,2))  
title('Offset around the jumps');
xlabel('Time around the jumps (sec)');
ylabel('Jump #');

%separate in two plots based on whether the jump was - or + 120 deg
jump_mag = shifted_offset(:,27);

figure('Position',[100 100 1000 800]),
subplot(1,2,1)
imagesc(shifted_offset(jump_mag<0,:))
colorbar
colormap(flipud(gray))
hold on
xline(26.5,'linewidth',2,'color','r');
xt = get(gca, 'XTick');                                            
xtlbl = linspace(min(short_time_AJ), max(short_time_AJ), numel(xt));                    
set(gca,'XTick',xt-2.5, 'XTickLabel',round(xtlbl,2))  
title('Offset around the jumps, negative changes');
xlabel('Time around the jumps (sec)');
ylabel('Jump #');

subplot(1,2,2)
imagesc(shifted_offset(jump_mag>0,:))
colorbar
colormap(flipud(gray))
hold on
xline(26.5,'linewidth',2,'color','r');
xt = get(gca, 'XTick');                                            
xtlbl = linspace(min(short_time_AJ), max(short_time_AJ), numel(xt));                    
set(gca,'XTick',xt-2.5, 'XTickLabel',round(xtlbl,2))  
title('Offset around the jumps, positive changes');
xlabel('Time around the jumps (sec)');
ylabel('Jump #');

saveas(gcf,[path,'\offset_AJ.png'])

%% sort by BM value pre cue jump

figure('Position',[100 100 1000 800]),
subplot(1,2,1)
imagesc(shifted_offset(I,:))
colorbar
gray_wrap = vertcat(gray,flipud(gray));
colormap(gray_wrap)
%colormap(flipud(gray))
hold on
xline(91,'linewidth',2,'color','r');
% xt = get(gca, 'XTick');                                            
% xtlbl = linspace(min(short_time_AJ), max(short_time_AJ), numel(xt));                    
% set(gca,'XTick',xt-2.5, 'XTickLabel',round(xtlbl,2))  
title('Offset around the jumps');
xlabel('Time around the jumps (sec)');
ylabel('Jump #');


% sort by BW value pre cue jump
subplot(1,2,2)
imagesc(shifted_offset(Iw,:))
colorbar
colormap(gray_wrap)
hold on
xline(91,'linewidth',2,'color','r');
% xt = get(gca, 'XTick');                                            
% xtlbl = linspace(min(short_time_AJ), max(short_time_AJ), numel(xt));                    
% set(gca,'XTick',xt-2.5, 'XTickLabel',round(xtlbl,2))  
title('Offset around the jumps');
xlabel('Time around the jumps (sec)');
ylabel('Jump #');

saveas(gcf,[path,'\offset_AJ_sorted_bm_and_bw.png'])

%% Look at shifted offset right pre-jump

shifted_offset_short = wrapTo180(offset-offset(:,90));
%short_time_AJ = time_around_jump(376:426);

figure('Position',[100 100 1000 800]),
subplot(1,2,1)
imagesc(shifted_offset_short(I,:))
colorbar
gray_wrap = vertcat(gray,flipud(gray));
colormap(gray_wrap)
hold on
xline(91,'linewidth',2,'color','r');
% xt = get(gca, 'XTick');                                            
% xtlbl = linspace(min(short_time_AJ), max(short_time_AJ), numel(xt));                    
% set(gca,'XTick',xt-2.5, 'XTickLabel',round(xtlbl,2))  
title({'Offset around the jumps','sorted by BM'});
xlabel('Time around the jumps');
ylabel('Jump #');


% sort by BW value pre cue jump
subplot(1,2,2)
imagesc(shifted_offset_short(Iw,:))
colorbar
colormap(gray_wrap)
hold on
xline(91,'linewidth',2,'color','r');
% xt = get(gca, 'XTick');                                            
% xtlbl = linspace(min(short_time_AJ), max(short_time_AJ), numel(xt));                    
% set(gca,'XTick',xt-2.5, 'XTickLabel',round(xtlbl,2))  
title({'Offset around the jumps','sorted by BW'});
xlabel('Time around the jumps');
ylabel('Jump #');

saveas(gcf,[path,'\shorter_offset_AJ_sorted_bm_and_bw.png'])

%% See how long it takes for the offset to go back to its pre-jump value

figure,
subplot(1,2,1)
plot(zscored_BM,first_frame,'o')
corr_value = corr(zscored_BM',first_frame, 'rows','complete');
hold on
text(0,1,['Corr = ',num2str(round(corr_value,2))]);
xlabel('Zscored bump magnitude');
ylabel('Zscored first frame post jump when offset returns');

subplot(1,2,2)
plot(zscored_BW,first_frame,'o')
corr_value = corr(zscored_BW',first_frame, 'rows','complete');
hold on
text(0,1,['Corr = ',num2str(round(corr_value,2))]);
xlabel('Zscored bump width');
ylabel('Zscored first frame post jump when offset returns');

saveas(gcf,[path,'\offset_full_returns.png'])


%% Look at this relationship when removing extreme values

first_frame_without_extremes = first_frame(first_frame<1);
zscored_BM_without_extremes = zscored_BM(first_frame<1);
zscored_BW_without_extremes = zscored_BW(first_frame<1);

figure,
subplot(1,2,1)
plot(zscored_BM_without_extremes,first_frame_without_extremes,'o')
corr_value2 = corr(zscored_BM_without_extremes',first_frame_without_extremes, 'rows','complete');
hold on
text(0,0,['Corr = ',num2str(corr_value2)]);
xlabel('Zscored bump magnitude');
ylabel('Zscored first frame post jump when offset returns');

subplot(1,2,2)
plot(zscored_BW_without_extremes,first_frame_without_extremes,'o')
corr_value2 = corr(zscored_BW_without_extremes',first_frame_without_extremes, 'rows','complete');
hold on
text(0,0,['Corr = ',num2str(corr_value2)]);
xlabel('Zscored bump width');
ylabel('Zscored first frame post jump when offset returns');

saveas(gcf,[path,'offset_full_return_without_extreme.png'])


%% Find nans (i.e., occasions when offset didn't return to previous status)

non_returns = find(isnan(first_frame));
%look at BM for those
figure,
plot(zscored_BM(non_returns),'ro')
ylabel('Zscored BM');
title('BM values for cue jumps where the offset never fully went back');
%it is very much all over the place.

saveas(gcf,[path,'\non_returns.png'])

%% Focus now on looking at the frames it takes to return 2/3 of the way...this implies looking at the first frame that's 40 deg away from the pre-jump offset

figure,
subplot(1,2,1)
plot(zscored_BM,partial_frame,'o')
corr_value_partial = corr(zscored_BM',partial_frame, 'rows','complete');
hold on
text(0,1,['Corr = ',num2str(corr_value_partial)]);
xlabel('Zscored bump magnitude');
ylabel('Zscored first frame post jump when offset returns 2/3');

subplot(1,2,2)
plot(zscored_BW,partial_frame,'o')
corr_value_partial = corr(zscored_BW',partial_frame, 'rows','complete');
hold on
text(0,1,['Corr = ',num2str(corr_value_partial)]);
xlabel('Zscored bump width');
ylabel('Zscored first frame post jump when offset returns 2/3');

saveas(gcf,[path,'partial_returns.png'])

%% Remove extreme values

partial_frame_without_extremes = partial_frame(partial_frame<1.5);
zscored_BM_without_extremes = zscored_BM(partial_frame<1.5);
zscored_BW_without_extremes = zscored_BW(partial_frame<1.5);

figure,
subplot(1,2,1)
plot(zscored_BM_without_extremes,partial_frame_without_extremes,'o')
corr_value_partial = corr(zscored_BM_without_extremes',partial_frame_without_extremes, 'rows','complete');
hold on
text(0,1,['Corr = ',num2str(corr_value_partial)]);
xlabel('Zscored bump magnitude');
ylabel('Zscored first frame post jump when offset returns 2/3');

subplot(1,2,2)
plot(zscored_BW_without_extremes,partial_frame_without_extremes,'o')
corr_value_partial = corr(zscored_BW_without_extremes',partial_frame_without_extremes, 'rows','complete');
hold on
text(0,1,['Corr = ',num2str(corr_value_partial)]);
xlabel('Zscored bump width');
ylabel('Zscored first frame post jump when offset returns 2/3');

saveas(gcf,[path,'partial_returns_without_extreme.png'])

%% Focus now on looking at the frames it takes to return 1/2 of the way...this implies looking at the first frame that's 60 deg away from the pre-jump offset

figure,
subplot(1,2,1)
plot(zscored_BM,half_frame,'o')
corr_value_half = corr(zscored_BM',half_frame, 'rows','complete');
hold on
text(0,1,['Corr = ',num2str(corr_value_half)]);
xlabel('Zscored bump magnitude');
ylabel('Zscored first frame post jump when offset returns 1/2');

subplot(1,2,2)
plot(zscored_BW,half_frame,'o')
corr_value_half = corr(zscored_BW',half_frame, 'rows','complete');
hold on
text(0,1,['Corr = ',num2str(corr_value_half)]);
xlabel('Zscored bump width');
ylabel('Zscored first frame post jump when offset returns 1/2');

saveas(gcf,[path,'half_returns.png'])


%% Remove extreme values

half_frame_without_extremes = half_frame(half_frame<1.5);
zscored_BM_without_extremes = zscored_BM(half_frame<1.5);
zscored_BW_without_extremes = zscored_BW(half_frame<1.5);

figure,
subplot(1,2,1)
plot(zscored_BM_without_extremes,half_frame_without_extremes,'o')
corr_value_partial = corr(zscored_BM_without_extremes',half_frame_without_extremes, 'rows','complete');
hold on
text(0,1,['Corr = ',num2str(corr_value_partial)]);
xlabel('Zscored bump magnitude');
ylabel('Zscored first frame post jump when offset returns 1/2');

subplot(1,2,2)
plot(zscored_BW_without_extremes,half_frame_without_extremes,'o')
corr_value_partial = corr(zscored_BW_without_extremes',half_frame_without_extremes, 'rows','complete');
hold on
text(0,1,['Corr = ',num2str(corr_value_partial)]);
xlabel('Zscored bump magnitude');
ylabel('Zscored first frame post jump when offset returns 1/2');

saveas(gcf,[path,'half_returns_without_extreme.png'])

%% Run model

allData = table(zscored_BM',zscored_BW',first_frame);
allData.Properties.VariableNames = {'zscored_BM','zscored_BW','first_frame'};

% Compute model for bump magnitude
mdl_BM = fitlme(allData,'first_frame~zscored_BM')

% Compute model for bump width
mdl_BW = fitlme(allData,'first_frame~zscored_BW')
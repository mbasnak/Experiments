%Code to generate a figure comparing accumulated heading and bump tracking
%in empty trials

clear all; close all;
folderNames = dir('Z:\Wilson Lab\Mel\Experiments\Uncertainty\Exp28\data');


%for the fly folders
for folder = 1:length(folderNames)
    if (contains(folderNames(folder).name,'60D05')==1 & contains(folderNames(folder).name,'20210127')==0 & contains(folderNames(folder).name,'20210218_60D05_7f_fly2')==0) %I'm excluding the first fly because I didn't run a darkness trial for it
        
        path = [folderNames(folder).folder,'\',folderNames(folder).name];
        %get the sessions info
        load([path,'\sessions_info.mat'])
        
        %load the empty trial and inverted gain trial data
        load([path,'\analysis\analysis_sid_',num2str(sessions_info.empty_trial),'_tid_0.mat']);
        
        nbins = 20;
        yaw_speed = abs(data.vel_yaw_ds);
        maxBinYS = min([max(yaw_speed),nanmean(yaw_speed)+3*nanstd(yaw_speed)]);
        binWidthYS = maxBinYS/nbins;
        YSBins = [0:binWidthYS:maxBinYS];
        for_vel = abs(data.vel_for_ds');
        maxBinFV = min([max(for_vel),nanmean(for_vel)+3*nanstd(for_vel)]);
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
        saveas(gcf,[path,'\analysis\plots\empty_trial_bm_heatmap.png']);
        
        bump_magnitude = data.bump_magnitude;
        save([path,'\analysis\vel_bm_data_et.mat'],'for_vel','yaw_speed','bump_magnitude');
        
    end
end
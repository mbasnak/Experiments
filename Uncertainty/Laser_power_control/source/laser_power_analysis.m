%Code to analyze the laser power control

%the idea behind this experiment is to see if differences in the laser
%power used to image the fly can result in differences in the bump
%magnitude and width at half max measured. We want to make sure that the
%differences in bump parameters we're seeing across conditions are not just
%experimental noise

%load the data
clear all; close all;

%list all of the files
path = 'Z:\Wilson Lab\Mel\Experiments\Uncertainty\Laser_power_control\data';
folderNames = dir(path);

for folder = 1:length(folderNames)
    if contains(folderNames(folder).name,'60D05')
        
        %load session info
        load(fullfile(folderNames(folder).folder,folderNames(folder).name,'sessions_info.mat'))
        
        %list the data file names
        fileNames = dir(fullfile(folderNames(folder).folder,folderNames(folder).name,'analysis'));
        
        %load the data corresponding to the different laser power values
        for file = 1:length(fileNames)
            if contains(fileNames(file).name,['sid_',num2str(sessions_info.laser_5)])
                data{folder,1} = load(fullfile(fileNames(file).folder,fileNames(file).name));
            elseif contains(fileNames(file).name,['sid_',num2str(sessions_info.laser_10)])
                data{folder,2} = load(fullfile(fileNames(file).folder,fileNames(file).name));
            elseif contains(fileNames(file).name,['sid_',num2str(sessions_info.laser_15)])
                data{folder,3} = load(fullfile(fileNames(file).folder,fileNames(file).name));
            end
        end
    end  
end

%remove empty cells
data = data(~cellfun(@isempty, data(:,1)), :);

%% Get bump parameters per laser power

figure,

%bump magnitude
for fly = 1:length(data)
    for laser_power = 1:3
        BM(fly,laser_power) = mean(data{fly,laser_power}.data.bump_magnitude);
    end
end

subplot(1,3,1)
plot(BM','-o','color',[.5 .5 .5])
hold on
errorbar(1:3,mean(BM),std(BM)/sqrt(length(BM)),'-ko');
xlim([0 4]);
xticks([1 2 3]);
xticklabels({'5','10','15'});
xlabel('Laser power (%)');
ylabel('Bump magnitude');


%bump width
for fly = 1:length(data)
    for laser_power = 1:3
        
        for timepoint = 1:size(data{fly,laser_power}.data.mean_dff_EB,2)
            
            %linearly interpolate to have 1000 datapoints instead of 8
            interp_ex_data = interp1([1:8],data{fly,laser_power}.data.mean_dff_EB(:,timepoint),[1:7/1000:8]);
            %Find the half max point
            half_max = (max(data{fly,laser_power}.data.mean_dff_EB(:,timepoint))-min(data{fly,laser_power}.data.mean_dff_EB(:,timepoint)))/2 + min(data{fly,laser_power}.data.mean_dff_EB(:,timepoint));
            [ex_bump_mag_interp I_interp] = max(interp_ex_data);
            %Find in each half the index closest to the half max
            diff_data = abs(interp_ex_data-half_max);
            [sortedVals,indexes] = sort(diff_data);
            %remove all the indexes that are less than 175 datapoints away (i.e., a little more than 1 glomerulus) from
            %index(1)
            diff_indexes = abs(indexes-indexes(1));
            indexes(diff_indexes<175 & diff_indexes>0)=NaN;
            indexes = indexes(~isnan(indexes));
            two_indexes = [indexes(1), indexes(2)];
            I1 = min(two_indexes);
            I2 = max(two_indexes);
            if (all(two_indexes>I_interp) | all(two_indexes<I_interp))
                half_max_w = I1+1000-I2;
            else
                half_max_w = I2-I1;
            end
            %convert to EB tiles
            half_max_width(timepoint) = half_max_w*8/1001;
        end
    BW(fly,laser_power) = mean(half_max_width);        
    end
end


subplot(1,3,2)
plot(BW','-o','color',[.5 .5 .5])
hold on
errorbar(1:3,mean(BW),std(BW)/sqrt(length(BW)),'-ko');
xlim([0 4]);
xticks([1 2 3]);
xticklabels({'5','10','15'});
xlabel('Laser power (%)');
ylabel('Bump width');

%total movement
for fly = 1:length(data)
    for laser_power = 1:3
        total_mvt(fly,laser_power) = nanmean(data{fly,laser_power}.data.total_mvt_ds);
    end
end

subplot(1,3,3)
plot(total_mvt','-o','color',[.5 .5 .5])
hold on
errorbar(1:3,mean(total_mvt),std(total_mvt)/sqrt(length(total_mvt)),'-ko');
xlim([0 4]);
xticks([1 2 3]);
xticklabels({'5','10','15'});
xlabel('Laser power (%)');
ylabel('Total movement (deg/s)');

saveas(gcf,[path,'\groupPlots\laser_control.png']);


%% Bin the data into different bump magnitude values and compute bump variability for each laser power, making one plot for each window size

%Define different bump magnitude bins
%1) Get all the bump magnitude data to find bin limits
allBumpMag = [];
for fly = 1:length(data)
    for laser_power = 1:size(data,2)
        allBumpMag = [allBumpMag,data{fly,laser_power}.data.bump_magnitude];
    end
end
%2) Set bin limits
minBin = prctile(allBumpMag,5);
maxBin = prctile(allBumpMag,95);
%3) Set bin width
nBins = 20;
bin_width = (maxBin-minBin)/nBins;
binLimits = [minBin:bin_width:maxBin];

%Define the different bin numbers
bins = [5,10,15,20,30,40,50,100];

for bin_size = 1:length(bins)
    
   figure,
   
   for laser_power = 1:size(data,2)
   %Bin bump magnitude values per laser power and get bump var per bin
        for bin = 1:nBins
            for fly = 1:length(data)
                bm = data{fly,laser_power}.data.bump_magnitude;
                bump_var_fly(fly,bin,laser_power) = circ_std(data{fly,laser_power}.data.dff_pva(bm > binLimits(bin) & bm < binLimits(bin+1))); 
            end
            bump_var = squeeze(nanmean(bump_var_fly));
        end   
   end
   
    plot(bump_var)   
    %xlim([minBin maxBin]);
    ylim([0 1.8]);
    
end


%% Bin the data to have the variability and bump magnitude computed over different
% bins

%define the possible bins
bins = [1,5,10,15,20,30,40,50,100];
clear all_binned_bm
clear all_binned_bv

%compute relationship between bump stability and bump parameters for each
%contrast and temporal window
for bin_amount = 1:length(bins)
    clear per_bin_bm 
    clear per_bin_bv
    clear bin_limits
    %figure,
    
    for fly = 1:length(data)
        for laser_power = 1:3
            
            %Get bin limits
            bin_number = bins(bin_amount);
            bin_width = floor(length(data{fly,laser_power}.data.dff_pva)/bin_number)-1;
            bin_limits{1} = [1,bin_width+1];
            for bin = 2:bin_number-1
                bin_limits{bin} = [2+bin_width*(bin-1),1+bin_width*(bin-1)+bin_width];
            end
            bin_limits{bin_number} = [length(data{fly,laser_power}.data.dff_pva)-bin_width,length(data{fly,laser_power}.data.dff_pva)];
            
            for bin = 1:bin_number
                %Get bump magnitude and bump variability per bin
                per_bin_bm{fly,laser_power,bin} = mean(data{fly,laser_power}.data.bump_magnitude(bin_limits{bin}(1):bin_limits{bin}(2)));
                [~, per_bin_bv{fly,laser_power,bin}] = circ_std(data{fly,laser_power}.data.dff_pva(bin_limits{bin}(1):bin_limits{bin}(2)));            
            end
                       
        end
            %Get mean binned bump magnitude and bump variability
            binned_bm = mean(cell2mat(per_bin_bm),3);
            binned_bv = mean(cell2mat(per_bin_bv),3);
            all_binned_bm{bin_amount} = mean(cell2mat(per_bin_bm),3);
            all_binned_bv{bin_amount} = mean(cell2mat(per_bin_bv),3);

    end
    
    %Get mean binned bump magnitude and bump variability across flies   
    %Plot
    figure,
    plot(mean(binned_bm),mean(binned_bv),'-ko')
    xlabel('Bump magnitude');
    ylabel('Bump variability');
    ylim([0.2 1.6]);
    title(['Bin number = ',num2str(bin_number)]);

end

%plot them all in one figure
figure('Position',[100 100 1000 800]),
for bin_amount = 1:length(bins)
    plot(mean(all_binned_bm{bin_amount}),mean(all_binned_bv{bin_amount}),'-o')
    hold on
end
legendCell = strcat('N bins=',string(num2cell(bins)));
legend(legendCell)
xlabel('Bump magnitude');
ylabel('Bump variability');
ylim([0.2 1.6]);



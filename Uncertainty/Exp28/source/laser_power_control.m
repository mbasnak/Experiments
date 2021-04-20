%bump stability and bump magnitude control

clear all; close all;

%% import data

%laser power 5%
lower = load('Z:\Wilson Lab\Mel\Experiments\Uncertainty\Exp28\data\20210127_60D05_7f\analysis\analysis_sid_3_tid_0.mat');
%laser power 10%
middle = load('Z:\Wilson Lab\Mel\Experiments\Uncertainty\Exp28\data\20210127_60D05_7f\analysis\analysis_sid_2_tid_0.mat');
%laser power 15%
upper = load('Z:\Wilson Lab\Mel\Experiments\Uncertainty\Exp28\data\20210127_60D05_7f\analysis\analysis_sid_4_tid_0.mat');

%% compute bump magnitude, width at half max, and bump stability for each dataset

%get bump position
data{1,1} = lower.data.dff_pva;
data{2,1} = middle.data.dff_pva;
data{3,1} = upper.data.dff_pva;

%compute bump magnitude
data{1,2} = max(lower.data.mean_dff_EB) - min(lower.data.mean_dff_EB);
data{2,2} = max(middle.data.mean_dff_EB) - min(middle.data.mean_dff_EB);
data{3,2} = max(upper.data.mean_dff_EB) - min(upper.data.mean_dff_EB); 


%define the possible bins
bins = [5,10,15,20,30,40,50,100];

%compute relationship between bump stability and bump parameters for each
%contrast and temporal window
for laser_power_level = 1:3
    for bin_amount = 1:length(bins)
        
        %reset variables
        clear binned_bump_mag
        clear corr_coefficient_bm
        clear pvalue_bm
        clear bump_var
        
        %change bin number
        bin_number = bins(bin_amount);
        bin_width = floor(length(data{laser_power_level,1})/bin_number);
        
        %get bin limits
        bin_limits{1} = [1,bin_width+1];
        for bin = 2:bin_number-1
            bin_limits{bin} = [2+bin_width*(bin-1),1+bin_width*(bin-1)+bin_width];
        end
        bin_limits{bin_number} = [length(data{laser_power_level,1})-bin_width,length(data{laser_power_level,1})];
        
        %get circular standard deviation of bump position
        for bin = 1:bin_number
            [~, bump_var(bin)] = circ_std(data{laser_power_level,1}(bin_limits{bin}(1):bin_limits{bin}(2)));
        end
        
        %get bump magnitude per bin
        for bin = 1:bin_number
            binned_bump_mag(bin) = median(data{laser_power_level,2}(bin_limits{bin}(1):bin_limits{bin}(2)));
        end
               
        %get correlation between variables
        [corr_coefficient_bm,pvalue_bm] = corrcoef(bump_var,binned_bump_mag);
        correlations_bm{laser_power_level}(bin_amount) = corr_coefficient_bm(1,2);
        pval_bm{laser_power_level}(bin_amount) = pvalue_bm(1,2);
        
        
    end
end

%Plot the relationships
figure('Position',[200 200 1200 800]),
for laser_power_level = 1:3
    subplot(3,1,laser_power_level)
    plot(bins(pval_bm{laser_power_level}<0.05),correlations_bm{laser_power_level}(pval_bm{laser_power_level}<0.05),'ro')
    hold on
    plot(bins(pval_bm{laser_power_level}>=0.05),correlations_bm{laser_power_level}(pval_bm{laser_power_level}>=0.05),'ko')
    line([0 100],[0 0],'color','b');
    ylim([-1 1]);
    xlabel('Bin number');
    ylabel('Correlation coefficient');
    title('Relationship between bump variation and bump magnitude');

end

%% get mean bump magnitude and bump variation

bump_mag = [mean(data{1,2}),mean(data{2,2}),mean(data{3,2})];
[~,bump_stb(1)] = circ_std(data{1,1});
[~,bump_stb(2)] = circ_std(data{2,1});
[~,bump_stb(3)] = circ_std(data{3,1});


figure,
plot(bump_mag,'-o')
hold on
plot(bump_stb,'-o')
legend('Bump magnitude','Bump instability')
xlim([0 4]);
xticks(1:3);
xticklabels({'5%','10%','15%'});
xlabel('Laser power');
%Bump stability analysis
clear all; close all;

%load example fly
load('Z:\Wilson Lab\Mel\Experiments\Uncertainty\Exp25\data\Experimental\two_ND_filters_3_contrasts\20201127_60D05_7f\analysis\analysis_sid_0_tid_0.mat')

%get bump velocity
bump_pos = data.dff_pva;

%Define bins
bin_number = 5;
bin_width = floor(length(bump_pos)/bin_number);
bin_limits{1} = [1,bin_width+1];
for bin = 2:bin_number-1
    bin_limits{bin} = [2+bin_width*(bin-1),1+bin_width*(bin-1)+bin_width];
end
bin_limits{bin_number} = [length(bump_pos)-bin_width,length(bump_pos)];


figure('Position', [200 100 1600 800]),
subplot(6,bin_number,[1 bin_number])
imagesc(data.dff_matrix)
colormap(gray)
hold on
for bin = 1:bin_number-1
    line([bin_width*bin bin_width*bin],[17 0],'color','r','lineWidth',2);
end
ylabel('PB glomerulus');
title('EPG activity');
set(gca,'xticklabel',{[]})

subplot(6,bin_number,[bin_number+1 bin_number*2])
plot(data.time,bump_pos)
hold on
for bin = 1:bin_number-1
    line([data.time(bin_width*bin) data.time(bin_width*bin)],[-pi pi],'color','r','lineWidth',2);
end
ylim([-pi pi]);
xlim([0 data.time(end)]);
title('Bump position');
ylabel('Radians');

%look at bump pos distribution in each section
for bin = 1:bin_number
    subplot(6,bin_number,bin_number*2+bin)
    polarhistogram(bump_pos(bin_limits{bin}(1):bin_limits{bin}(2))) 
    Ax = gca; % current axes
    Ax.ThetaGrid = 'off';
    Ax.RGrid = 'off';
    Ax.RTickLabel = [];
    Ax.ThetaTickLabel = [];
end

%get circular standard deviation of bump position
for bin = 1:bin_number
    [~, bump_var(bin)] = circ_std(bump_pos(bin_limits{bin}(1):bin_limits{bin}(2)));
end
%plot
subplot(6,bin_number,[bin_number*3+1 bin_number*4])
plot(bump_var,'-o');
xlim([0.5 bin_number+0.5]);
title('Bump position variation');
ylabel('Circ_std');
xticks([1 2 3 4 5])
set(gca,'xticklabel',{[]})

% Get bump magnitude and relate to bump stability
%compute bump magnitude as max-min
bump_mag = max(data.mean_dff_EB)-min(data.mean_dff_EB); 
%get bump magnitude per bin
for bin = 1:bin_number
    binned_bump_mag(bin) = mean(bump_mag(bin_limits{bin}(1):bin_limits{bin}(2)));
end
subplot(6,bin_number,[bin_number*4+1 bin_number*5])
plot(binned_bump_mag,'-o');
xlim([0.5 bin_number+0.5]);
title('Mean bump magnitude');
ylabel('Bump mag (max-min)');
xticks([1 2 3 4 5])
set(gca,'xticklabel',{[]})

% Get bump width at half max and relate to bump stability
%compute bump width at half max
for timepoint = 1:size(data.mean_dff_EB,2)
    
    %linearly interpolate to have 1000 datapoints instead of 8
    interp_ex_data = interp1([1:8],data.mean_dff_EB(:,timepoint),[1:7/1000:8]);
    %Find the half max point
    half_max = (max(data.mean_dff_EB(:,timepoint))-min(data.mean_dff_EB(:,timepoint)))/2 + min(data.mean_dff_EB(:,timepoint));
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

%get bump width at half max per bin
for bin = 1:bin_number
    binned_half_width(bin) = median(half_max_width(bin_limits{bin}(1):bin_limits{bin}(2)));
end
subplot(6,bin_number,[bin_number*5+1 bin_number*6])
plot(binned_half_width,'-o');
xlim([0.5 bin_number+0.5]);
title('Mean bump width at half max');
ylabel('BWHM');
xticks([1 2 3 4 5])
set(gca,'xticklabel',{[]})

%% Get elationship between bump parameters and bump position variation

figure,
subplot(1,2,1)
scatter(bump_var,binned_bump_mag)
xlabel('Bump position variation (circ_std)');
ylabel('Bump magnitude (max-min)');
%get correlation between variables
[corr_coefficient,pvalue] = corrcoef(bump_var,binned_bump_mag);
corr_coef_1 = corr_coefficient(1,2);
pval_1 = pvalue(1,2);
text(max(bump_var)-0.6,max(binned_bump_mag),['r = ',num2str(round(corr_coef_1,2))]);
text(max(bump_var)-0.6,max(binned_bump_mag)-0.05,['p-val = ',num2str(round(pval_1,2))]);

subplot(1,2,2)
scatter(bump_var,binned_half_width)
xlabel('Bump position variation (circ_std)');
ylabel('Bump width at half max');
%get correlation between variables
[corr_coefficient,pvalue] = corrcoef(bump_var,binned_half_width);
corr_coef_2 = corr_coefficient(1,2);
pval_2 = pvalue(1,2);
text(max(bump_var)-0.6,max(binned_half_width),['r = ',num2str(round(corr_coef_2,2))]);
text(max(bump_var)-0.6,max(binned_half_width)-0.05,['p-val = ',num2str(round(pval_2,2))]);



%% Vary the bin number

clear binned_bump_mag
clear binned_half_width
clear corr_coefficient_bm
clear corr_coefficient_hw
clear pvalue_bm
clear pvalue_hw
clear bump_var

bins = [5,10,15,20,30,40,50,100];

for bin_amount = 1:length(bins)
    
    clear binned_bump_mag
    clear binned_half_width
    clear corr_coefficient_bm
    clear corr_coefficient_hw
    clear pvalue_bm
    clear pvalue_hw
    clear bump_var
    
    %change bin number
    bin_number = bins(bin_amount);
    bin_width = floor(length(bump_pos)/bin_number);
    
    %get bin limits
    bin_limits{1} = [1,bin_width+1];
    for bin = 2:bin_number-1
        bin_limits{bin} = [2+bin_width*(bin-1),1+bin_width*(bin-1)+bin_width];
    end
    bin_limits{bin_number} = [length(bump_pos)-bin_width,length(bump_pos)];
    
    %get circular standard deviation of bump position
    for bin = 1:bin_number
        [~, bump_var(bin)] = circ_std(bump_pos(bin_limits{bin}(1):bin_limits{bin}(2)));
    end
    
    %get bump magnitude per bin
    for bin = 1:bin_number
        binned_bump_mag(bin) = median(bump_mag(bin_limits{bin}(1):bin_limits{bin}(2)));
    end
    
    %get bump width at half max per bin
    for bin = 1:bin_number
        binned_half_width(bin) = median(half_max_width(bin_limits{bin}(1):bin_limits{bin}(2)));
    end
    
    %get correlation between variables
    [corr_coefficient_bm,pvalue_bm] = corrcoef(bump_var,binned_bump_mag);
    corr_coef_bm(bin_amount) = corr_coefficient_bm(1,2);
    pval_bm(bin_amount) = pvalue_bm(1,2);    
    
    [corr_coefficient_hw,pvalue_hw] = corrcoef(bump_var,binned_half_width);
    corr_coef_hw(bin_amount) = corr_coefficient_hw(1,2);
    pval_hw(bin_amount) = pvalue_hw(1,2);

end

%Plot the relationships
figure('Position',[200 200 1200 800]),
subplot(1,2,1)
plot(bins(pval_bm<0.05),corr_coef_bm(pval_bm<0.05),'ro')
hold on
plot(bins(pval_bm>=0.05),corr_coef_bm(pval_bm>=0.05),'ko')
line([1 100],[0 0],'color','b');
ylim([-1 1]);
xticks([5,10,15,20,30,40,50,100]);
xticklabels({'5','10','15','20','30','40','50','100'});
%Add custom legend
h = zeros(2, 1);
h(1) = plot(NaN,NaN,'ro');
h(2) = plot(NaN,NaN,'ko');
legend(h, 'Significant','Not significant');
xlabel('Bin number');
ylabel('Correlation coefficient');
title('Relationship between bump variation and bump magnitude');

subplot(1,2,2)
plot(bins(pval_hw<0.05),corr_coef_hw(pval_hw<0.05),'ro')
hold on
plot(bins(pval_hw>=0.05),corr_coef_hw(pval_hw>=0.05),'ko')
line([1 100],[0 0],'color','b');
xlabel('Bin number');
ylim([-1 1]);
%Add custom legend
h = zeros(2, 1);
h(1) = plot(NaN,NaN,'ro');
h(2) = plot(NaN,NaN,'ko');
legend(h, 'Significant','Not significant');
xticks([5,10,15,20,30,40,50,100]);
xticklabels({'5','10','15','20','30','40','50','100'});
ylabel('Correlation coefficient');
title('Relationship between bump variation and bump width at half max');

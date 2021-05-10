new_order = [16,7,15,6,14,5,13,4,12,3,11,2,10,1,9,8];

figure('Position',[100 100 1000 800]),
subplot(3,1,1)
plot(data.dff_matrix(:,1000),'-o')
xlim([1 16]);
xticks([1:16])
xticklabels({'G8R','G7R','G6R','G5R','G4R','G3R','G2R','G1R','G1L','G2L','G3L','G4L','G5L','G6L','G7L','G8L'})
ylabel('DF/F');
title('Data from one end of the PB to the other');

subplot(3,1,2)
plot(data.shifted_dff_matrix(:,1000),'-o')
xlim([1 16]);
xticks([1:16])
xticklabels({'G1L','G8L','G7L','G6L','G5L','G4L','G3L','G2L','G1R','G2R','G3R','G4R','G5R','G6R','G7R','G8R'})
ylabel('DF/F');
title('Shifted data such that both halves can be averaged (G1 moves)');

subplot(3,1,3)
plot(data.dff_matrix(new_order,1000),'-o')
xlim([1 16]);
xticks([1:16])
xticklabels({'G8L','G2R','G7L','G3R','G6L','G4R','G5L','G5R','G4L','G6R','G3L','G7R','G2L','G8R','G1L','G1R'})
ylabel('DF/F');
title('Data sorted to match EB coordinates');

saveas(gcf,['Z:\Wilson Lab\Mel\Experiments\Uncertainty\Exp28\data\20210127_60D05_7f\analysis\plots\sort_dff_data_ex3.png']);


%% Compute half width with the new coordinates

test_data = data.dff_matrix(new_order,:);

for timepoint = [1,100,1000]
    
    ex_bump_mag = max(test_data(:,timepoint)) - min(test_data(:,timepoint));
    %Linearly interpolate the data in 'EB coordinates' to have 1000 datapoints instead of 8
    interp_ex_data = interp1([1:16],test_data(:,timepoint),[1:15/1000:16]);
    %Find the half max value in the y axis as the middle point between the min
    %and max y axis values
    half_max = (max(test_data(:,timepoint))-min(test_data(:,timepoint)))/2 + min(test_data(:,timepoint));
    [ex_bump_mag_interp I_interp] = max(interp_ex_data);
    %Find in each half the index closest to the half max value
    diff_data = abs(interp_ex_data-half_max);
    [sortedVals,indexes] = sort(diff_data);
    %Remove all the indexes that are less than 175 datapoints away (i.e., a little more than 1 glomerulus) from
    %index(1)
    diff_indexes = abs(indexes-indexes(1));
    indexes(diff_indexes<175 & diff_indexes>0)=NaN;
    indexes = indexes(~isnan(indexes));
    two_indexes = [indexes(1), indexes(2)];
    I1 = min(two_indexes);
    I2 = max(two_indexes);
    %Convert to EB coordinates
    half_max_width(timepoint) = half_max_w*16/1001;
    
    figure,
    plot(interp_ex_data)
    hold on
    plot([I_interp I_interp],[ex_bump_mag_interp ex_bump_mag_interp],'ko')
    plot(I1,interp_ex_data(I1),'ro')
    plot(I2,interp_ex_data(I2),'ro')
    if (all(two_indexes>I_interp) | all(two_indexes<I_interp))
        line([1 I1],[half_max half_max],'color','r','LineWidth',2);
        line([I2 1000],[half_max half_max],'color','r','LineWidth',2);
        half_max_w = I1+1000-I2;
    else
        line([I1 I2],[half_max half_max],'color','r','LineWidth',2);
        half_max_w = I2-I1;
    end
    text(500,ex_bump_mag,num2str(half_max_w))
    
    saveas(gcf,['Z:\Wilson Lab\Mel\Experiments\Uncertainty\Exp28\data\20210127_60D05_7f\analysis\plots\half_width_ex',num2str(timepoint),'.png']);

end


%% Compute half width with the usual coordinates

test_data2 = data.mean_dff_EB(:,:);

for timepoint = [1,100,1000]
    
    ex_bump_mag = max(test_data2(:,timepoint)) - min(test_data2(:,timepoint));
    %Linearly interpolate the data in 'EB coordinates' to have 1000 datapoints instead of 8
    interp_ex_data = interp1([1:8],test_data2(:,timepoint),[1:7/1000:8]);
    %Find the half max value in the y axis as the middle point between the min
    %and max y axis values
    half_max = (max(test_data2(:,timepoint))-min(test_data2(:,timepoint)))/2 + min(test_data2(:,timepoint));
    [ex_bump_mag_interp I_interp] = max(interp_ex_data);
    %Find in each half the index closest to the half max value
    diff_data = abs(interp_ex_data-half_max);
    [sortedVals,indexes] = sort(diff_data);
    %Remove all the indexes that are less than 175 datapoints away (i.e., a little more than 1 glomerulus) from
    %index(1)
    diff_indexes = abs(indexes-indexes(1));
    indexes(diff_indexes<175 & diff_indexes>0)=NaN;
    indexes = indexes(~isnan(indexes));
    two_indexes = [indexes(1), indexes(2)];
    I1 = min(two_indexes);
    I2 = max(two_indexes);
    %Convert to EB coordinates
    half_max_width(timepoint) = half_max_w*8/1001;
    
    figure,
    plot(interp_ex_data)
    hold on
    plot([I_interp I_interp],[ex_bump_mag_interp ex_bump_mag_interp],'ko')
    plot(I1,interp_ex_data(I1),'ro')
    plot(I2,interp_ex_data(I2),'ro')
    if (all(two_indexes>I_interp) | all(two_indexes<I_interp))
        line([1 I1],[half_max half_max],'color','r','LineWidth',2);
        line([I2 1000],[half_max half_max],'color','r','LineWidth',2);
        half_max_w = I1+1000-I2;
    else
        line([I1 I2],[half_max half_max],'color','r','LineWidth',2);
        half_max_w = I2-I1;
    end
    text(500,ex_bump_mag,num2str(half_max_w))
    
    saveas(gcf,['Z:\Wilson Lab\Mel\Experiments\Uncertainty\Exp28\data\20210127_60D05_7f\analysis\plots\usual_half_width_ex',num2str(timepoint),'.png']);

    
end
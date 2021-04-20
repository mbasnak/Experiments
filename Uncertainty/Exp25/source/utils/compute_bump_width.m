%function to compute the width at half max of the bump at a specific
%timepoint

%OUTPUT = with at half max value at a specific timepoint
%INPUTS
    %data = dff data in EB coordinates
    %timepoint = point in time at which we're evaluating the half width

function [half_max_width] = compute_bump_width(data,timepoint)
    
    %Linearly interpolate the data in 'EB coordinates' to have 1000 datapoints instead of 8
    interp_ex_data = interp1([1:8],data(:,timepoint),[1:7/1000:8]);
    %Find the half max value in the y axis as the middle point between the min
    %and max y axis values
    half_max = (max(data(:,timepoint))-min(data(:,timepoint)))/2 + min(data(:,timepoint));
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
    if (all(two_indexes>I_interp) | all(two_indexes<I_interp))
        half_max_w = I1+1000-I2;  
    else
        half_max_w = I2-I1;
    end  
    %Convert to EB coordinates
    half_max_width = half_max_w*8/1001;     
end
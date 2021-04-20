for i = 1:num_rois %for each ROI
    %% Look for the rois that match the roi num
    indices = find(strcmp({roi(:).name},char(roi_names(i)))); %get indices for the different ROIs that match a name
    summed = [];
    pixels = zeros(1, length(indices));
    
    %% Sum the values across all the rois
    for j = 1:length(indices)
        index = indices(j);
        summedData = squeeze(sum(regProductGCamp,3));
        roi_im = bsxfun(@times, squeeze(summedData(:,:,:)), roi(index).BW);
        summed = [summed squeeze(sum(sum(roi_im,1)))];
        mask = roi(index).BW; %mask of positions for the ROI in question
        pixels(j) = sum(mask(:)); %number of pixels that comprise that ROI
        
    end
    %% Calculate the trace as the mean value in the rois (by pixel number)
    figure,
    suptitle(['ROI #',num2str(i)]);
    trace = summed*pixels'./sum(pixels);
    subplot 321, plot(trace);
    [sorted,idx] = sort(trace); %sort the mean fluorescence values from lowest to highest
    subplot 322, plot(sorted);
    timepoints = size(trace,1); %get the number of volumes (which will be timepoints)
    
    %% Detect if "frames lost"
    
    test = diff(sorted);
    subplot 323, plot(test)
    values = find(test(1:end-100)>0.18E4); %this value sometimes needs to be tweeked.
    if (i~=4)
        if length(values)>0
            inflexionPoint = trace(idx(values(1)+1));
            trace(trace<inflexionPoint) = nan;
            sorted = sort(trace(~isnan(trace)));
        else
        end
    end

    subplot 324, plot(trace)
    subplot 325, plot(sorted)
    
    %% Baseline is the bottom fifth percentile
    
    fifthpercentile = floor(.05*timepoints);
    baseline_f = mean(sorted(1:fifthpercentile));
    
    %% Calculate dff and zscore
    dff = (trace-baseline_f)./baseline_f;
    dff = fillmissing(dff,'linear');
    subplot 326, plot(dff)
    trace = fillmissing(trace,'linear');
    z = zscore(trace);
    
end
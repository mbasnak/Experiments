%code to compute several ROI characteristics, among them dff/f and z-score.
%made by JL

function roi_data = ROI_analysis(parentDir, sid, tid_list)

%% ROI Analysis Function
% Objectives: From ROIs given, create an ROI data structure with following
% information:
%
% roi_analysis.baseline_f
% roi_analysis.dff
% roi_analysis.zscore
% roi_analysis.num
% roi_analysis.xi
% roi_analysis.yi
% roi_analysis.zplane
% roi_analysis.leftPB
% roi_analysis.rightPB
% roi_analysis.glomerulus
% roi_analysis.pixels
%

%% Look for session Files post registration. Find if there is one or two
% channels.
global slash;
if isunix() == 1
    slash = '/';
else
    slash = '\';
end
imagingDir = [parentDir slash '2p' slash];
num_tids = length(tid_list);

for i = 1:num_tids %for every trial
    
    %% Get files
    tid = tid_list(i); %get the trial id
    imagingTrialDir = [imagingDir 'sid_' num2str(sid) '_tid_' num2str(tid) slash]; %set the directory where the registered images live
    if ~isdir(imagingTrialDir)
        disp('Dir does not exist');
        return;
    end
    expression = ['*sid_' num2str(sid) '_tid_' num2str(tid) '*'];
    expression2 = ['rigid_*'];
    roiDir = [imagingDir slash 'ROI' slash]; %get the directory where the ROIs are saved
    expression3 = ['*sid_' num2str(sid) '_tid_' num2str(tid) '*']; %this looks just like 'expression'
    imagingFile = dir(fullfile(imagingTrialDir, expression2)); %get the rigid registration file info (all the data) from the corresponding trial folder
    roiFile = dir(fullfile(roiDir, expression3)); %get the ROI file info for the current trial from the ROI folder

    numChannels = size(imagingFile, 1);

    %% Load the GCaMP stack of each channel.
    
    summedStack = [];
    numSlices = [];
    numVols = [];

    if numChannels == 1
        tic;
        load(fullfile(imagingTrialDir, imagingFile.name)); %load the registered stack
        numSlices = size(regProduct, 3); %get the number of slices per volume (it'll usually be 8 for me: the original 12 - 4 flyback ones)
        numVols = size(regProduct, 4); %get the number of volumes (total times the numSlices have been imaged)
        summedStack = max(regProduct,[], 3); %this gets the maximum fluorescence value of the 8 slices for each volume imaged
        regProductGCamp = regProduct; %we're only using 1 channel so the GCaMP data is our channel 1 data
        toc;
    else % numChannels == 2
        load(fullfile(imagingTrialDir, imagingFile(1).name));
        summedStack = max(regProduct,[], 3);
        numSlices = size(regProduct, 3);
        numVols = size(regProduct, 4);

        regProductGCamp = regProduct;
        tic;
        load(fullfile(imagingTrialDir, imagingFile(2).name));
        regProductTdTom = regProduct;
        toc;
    end

    tic;
    %% Load the ROI file
    load(fullfile(roiDir, roiFile(1).name)); % loads variable roi, that has the position info about the ROIs

    toc;
    %% For each ROI
    roi_names = unique({roi(:).name}); %returns the names of the different ROIs, without repeting them if some are repeated
    num_rois = length(roi_names); %how many ROIs you have

    roi_data = [];

    tic;
    for i = 1:num_rois %for each ROI
       %% Look for the rois that match the roi num
       indices = find(strcmp({roi(:).name},char(roi_names(i)))); %get indices for the different ROIs that match a name
       summed = [];
       pixels = zeros(1, length(indices));
       
       %% Sum the values across all the rois
       for j = 1:length(indices)
           index = indices(j);
           
           %roi_im = bsxfun(@times, squeeze(regProductGCamp(:,:,roi(index).z,:)), roi(index).BW);
           %the function 'times' multiples matrices element-wise.
           %using bsxfun applied the function, in this case 'times' to the
           %matrices mentioned later.
           %we are then multiplying the GCaMP image (in this case slice 1)
           %by the positions of the current ROI we're assessing, such that
           %we obtain fluorescence values for those pixels we're interested
           %in
           
           %I'm trying to use the summed ROI activity instead of just the
           %max z slice so I commented the previous roi_im and added the following - MB 20200219
           summedData = squeeze(sum(regProductGCamp,3));
           roi_im = bsxfun(@times, squeeze(summedData(:,:,:)), roi(index).BW);
           
           summed = [summed squeeze(sum(sum(roi_im,1)))];
           %gives you back for each ROI a vector with a number of elements
           %equal to the number of volumes, where each element has the sum
           %of the total ROI fluorescence for that volume (for z 1 in this
           %case)
           
           mask = roi(index).BW; %mask of positions for the ROI in question
           pixels(j) = sum(mask(:)); %number of pixels that comprise that ROI
           %pixels seem to end up with one dimension. is the j in case
           %there are several z slices we're considering?
       end
       %% Calculate the trace as the mean value in the rois (by pixel number)
       trace = summed*pixels'./sum(pixels); %I don't really get this, trace is the same as summed
       %shouldn't this just be divided by the pixels, withour the
       %multiplication?
       %I think this makes sense if we're using multiple z slices
       
       [sorted,idx] = sort(trace); %sort the mean fluorescence values from lowest to highest
       timepoints = size(trace,1); %get the number of volumes (which will be timepoints)
       
       %% Detect if "frames lost"
  
       %I made a new code to detect lost frames - MB 20200128
       test = diff(sorted);
       values = find(test(1:end-100)>0.35E4);
       if length(values)>0
            inflexionPoint = trace(idx(values(1)+1));
            trace(trace<inflexionPoint) = nan;
            sorted = sort(trace(~isnan(trace)));
       else
       end
       
      %how is this deciding that the frames are lost? what is this actually
      %clustering? is it 'high' and 'low' response? and too low might mean
      %lost?
%        if sorted(end)/sorted(1) >= 30 %why is this criterion taken?
%            T = clusterdata(trace, 2); %this clusters the mean data in 2 groups
%            if mode(T) == 1 %if there are more observations from cluster 1
%                bad_cluster = 2; %the clusterdata 2 is 'bad'
%            else %if there are more observations from cluster2
%                bad_cluster = 1; %the cluster 1 is 'bad'
%            end
%            trace(T == bad_cluster) = nan; %this step is adding a lot of NaNs to 'lost' frames
%            sorted = sort(trace(T == 2)); 
%            timepoints = size(sorted, 1); %obtain the new timepoints (taking out the 'lost' frames
%        end
       
%Alternative code por removing lost frames - MB 20200123

%        if sorted(end)/sorted(1) >= 30 %why is this criterion taken?
%            T = clusterdata(trace, 2); %this clusters the mean data in 2 groups
%            if (median(trace(T==1)) > median(trace(T==2))) %& (sum(T==1) > 0.9*timepoints) 
%                bad_cluster = 2; %the clusterdata 2 is 'bad'
%            else
%                bad_cluster = 1;
%            end
%            trace(T == bad_cluster) = nan; %this step is adding a lot of NaNs to 'lost' frames
%            sorted = sort(trace(T ~= bad_cluster)); 
%            timepoints = size(sorted, 1); %obtain the new timepoints (taking out the 'lost' frames
%        end



       %% Baseline is the bottom fifth percentile
       %this is taking the baseline activity as the mean of the first 5% values of
       %timepoints, and since we pre-sorted them, that corresponds to the
       %lowest 5% values of the trial
       
       fifthpercentile = floor(.05*timepoints);
       baseline_f = mean(sorted(1:fifthpercentile));
       
       %% Calculate dff and zscore
       dff = (trace-baseline_f)./baseline_f;     
       
       %I'm adding this to linearly interpolate the missing values - MB
       %20200123
       dff = fillmissing(dff,'linear');
       trace = fillmissing(trace,'linear');
       z = zscore(trace);
       %z = (trace-nanmean(trace))./nanstd(trace);
       
       %% Create roi_analysis object for this ROI
       
       %save the analysis data obtained in a struct called roi_data
       roi_analysis = [];
       roi_analysis.baseline_f = baseline_f;
       roi_analysis.dff = dff;
       roi_analysis.zscore = z;
       %roi_analysis.num = i;
       roi_analysis.xi = {roi(indices).xi};
       roi_analysis.yi = {roi(indices).yi};
       roi_analysis.zplane = [roi(indices).z];
       %roi_analysis.leftPB = roi(indices(1)).leftPB;
       %roi_analysis.rightPB = roi(indices(1)).rightPB;
       %roi_analysis.glomerulus = roi(indices(1)).glomerulus;
       roi_analysis.name = char(roi_names(i));
       mask = roi(i).BW;
       roi_analysis.pixels = sum(mask(:));

       %% 
       if numChannels == 2
           summed_tdTom = [];
           pixels_tdTom = zeros(1, length(indices));
           for j = 1:length(indices)
               index = indices(j);
               roi_im = bsxfun(@times, squeeze(regProductTdTom(:,:,roi(index).z,:)), roi(index).BW);
               summed_tdTom = [summed_tdTom squeeze(sum(sum(roi_im,1)))];
               mask = roi(index).BW;
               pixels_tdTom(j) = sum(mask(:));
           end
           trace_tdTom = summed_tdTom*pixels_tdTom'./sum(pixels_tdTom);
           sorted_tdTom = sort(trace_tdTom);
           timepoints_tdTom = size(trace_tdTom,1);
           fifthpercentile_tdTom = floor(.05*timepoints_tdTom);
           baseline_f_tdTom = mean(sorted_tdTom(1:fifthpercentile_tdTom));
           dff_tdTom = (trace_tdTom-baseline_f_tdTom)./baseline_f_tdTom;
           z_tdTom = zscore(trace_tdTom);
           roi_analysis.baseline_f_tdTom = baseline_f_tdTom;
           roi_analysis.dff_tdTom = dff_tdTom;
           roi_analysis.zscore_tdTom = z_tdTom;
       end
       roi_data = [roi_data roi_analysis];
    end

    toc;
    tic;

    % Save data
    roi_analysis_path = [parentDir slash '2p' slash 'ROI_analysis'];
    if(~exist(roi_analysis_path, 'dir'))
    mkdir(roi_analysis_path); %generates a new directory called ROI_analysis
    end
    filename = [roi_analysis_path slash ['ROI_data_sid_' num2str(sid) '_tid_' num2str(tid) '.mat']];
    
    save(filename, 'roi_data'); %saves the file with the ROI analysis for the current trial
    toc;
end

%disp('ROI analysis done!');
end
%code to compute several the DF/F in several different ways

function dff_data = DFF_methods(parentDir, sid, tid_list)



%% Look for session Files post registration.

imagingDir = [parentDir '\2p\'];
num_tids = length(tid_list);

for i = 1:num_tids %for every trial
    
    %% Get files
    tid = tid_list(i); %get the trial id
    imagingTrialDir = [imagingDir 'sid_' num2str(sid) '_tid_' num2str(tid) '\']; %set the directory where the registered images live
    if ~isdir(imagingTrialDir)
        disp('Dir does not exist');
        return;
    end
    expression = ['*sid_' num2str(sid) '_tid_' num2str(tid) '*'];
    expression2 = ['rigid_*'];
    roiDir = [imagingDir '\ROI\']; %get the directory where the ROIs are saved
    expression3 = ['*sid_' num2str(sid) '_tid_' num2str(tid) '*']; %this looks just like 'expression'
    imagingFile = dir(fullfile(imagingTrialDir, expression2)); %get the rigid registration file info (all the data) from the corresponding trial folder
    roiFile = dir(fullfile(roiDir, expression3)); %get the ROI file info for the current trial from the ROI folder

    numChannels = size(imagingFile, 1);

    %% Load the GCaMP stack of each channel.
    
    summedStack = [];
    numSlices = [];
    numVols = [];
   
    load(fullfile(imagingTrialDir, imagingFile.name)); %load the registered stack
    numSlices = size(regProduct, 3); %get the number of slices per volume (it'll usually be 8 for me: the original 12 - 4 flyback ones)
    numVols = size(regProduct, 4); %get the number of volumes (total times the numSlices have been imaged)
    summedStack = max(regProduct,[], 3); %this gets the maximum fluorescence value of the 8 slices for each volume imaged
    regProductGCamp = regProduct; %we're only using 1 channel so the GCaMP data is our channel 1 data

    %% Load the ROI file
    load(fullfile(roiDir, roiFile(1).name)); % loads variable roi, that has the position info about the ROIs

    %% For each ROI
    roi_names = unique({roi(:).name}); %returns the names of the different ROIs, without repeting them if some are repeated
    num_rois = length(roi_names); %how many ROIs you have

    roi_data = [];


    for i = 1:num_rois %for each ROI
       %% Look for the rois that match the roi num
       indices = find(strcmp({roi(:).name},char(roi_names(i)))); %get indices for the different ROIs that match a name
       summed = [];
       pixels = zeros(1, length(indices));
       
       %% Sum the values across all the rois
       for j = 1:length(indices)
           index = indices(j); 
           
           %Get the summed ROI activity
           summedData = squeeze(sum(regProductGCamp,3));
           roi_im = bsxfun(@times, squeeze(summedData(:,:,:)), roi(index).BW);
           summed = [summed squeeze(sum(sum(roi_im,1)))];
           %gives you back for each ROI a vector with a number of elements
           %equal to the number of volumes, where each element has the sum
           %of the total ROI fluorescence for that volume (for z 1 in this
           %case)
           mask = roi(index).BW; %mask of positions for the ROI in question
           pixels(j) = sum(mask(:)); %number of pixels that comprise that ROI
           %pixels seem to end up with one dimension.
           
           
           %Get the average ROI activity
           
           
           
       end
       %% Calculate the trace as the mean value in the rois (by pixel number)
       trace = summed*pixels'./sum(pixels); %This makes sense if we're using multiple z slices
       [sorted,idx] = sort(trace); %sort the mean fluorescence values from lowest to highest
       timepoints = size(trace,1); %get the number of volumes (which will be timepoints)
       
       %% Detect if "frames lost"
  
       %I made a new code to detect lost frames - MB 20200128
       test = diff(sorted);
       values = find(test(1:end-100)>0.18E4); %this value sometimes needs to be tweeked.
   
       if length(values)>0
           inflexionPoint = trace(idx(values(1)+1));
           trace(trace<inflexionPoint) = nan;
           sorted = sort(trace(~isnan(trace)));
       else
       end
       
       %Determine percentage of frames lost
       percent_lost = sum(trace == nan)/length(trace);
       
       %% Compute the baseline
       
       %We are going to use different methods to compute the baseline
       
       %Fifth percentile
       fifthpercentile = floor(.05*timepoints);
       baseline_f_5 = mean(sorted(1:fifthpercentile));
       
       %Tenth percentile
       tenthpercentile = floor(.1*timepoints);
       baseline_f_10 = mean(sorted(1:tenthpercentile));     
       
       
       
       %% Calculate dff and zscore
       
       dff_5 = (trace-baseline_f_5)./baseline_f_5; 
       dff_10 = (trace-baseline_f_10)./baseline_f_10;        
       
       %Linearly interpolate the missing values
       dff_5 = fillmissing(dff_5,'linear');
       dff_10 = fillmissing(dff_10,'linear');
       trace = fillmissing(trace,'linear');
       
       %Get the z-scored data
       z = zscore(trace);
       
       %% Create roi_analysis object for this ROI
       
       %save the analysis data obtained in a struct called roi_data
       roi_analysis = [];
       roi_analysis.baseline_f = baseline_f;
       roi_analysis.dff = dff;
       roi_analysis.zscore = z;
       roi_analysis.trace = trace;
       roi_analysis.xi = {roi(indices).xi};
       roi_analysis.yi = {roi(indices).yi};
       roi_analysis.zplane = [roi(indices).z];
       roi_analysis.name = char(roi_names(i));
       mask = roi(i).BW;
       roi_analysis.pixels = sum(mask(:));


       roi_data = [roi_data roi_analysis];
    end


    % Save data
%     roi_analysis_path = [parentDir slash '2p' slash 'ROI_analysis'];
%     if(~exist(roi_analysis_path, 'dir'))
%     mkdir(roi_analysis_path); %generates a new directory called ROI_analysis
%     end
%     filename = [roi_analysis_path slash ['ROI_data_sid_' num2str(sid) '_tid_' num2str(tid) '.mat']];
%     
%     save(filename, 'roi_data'); %saves the file with the ROI analysis for the current trial
end


end
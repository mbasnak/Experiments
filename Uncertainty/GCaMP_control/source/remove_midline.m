function remove_midline(path)

%Code to remove the midline from the ROI selection and save as the ROIs for
%the old analysis

ROI_files = dir([path,'\2p\ROI']);

for file = 1:length(ROI_files)
    if contains(ROI_files(file).name,'midline')
        load([path,'\2p\ROI\',ROI_files(file).name])
        roi(17) = [];
        save([path,'\2p\ROI\',ROI_files(file).name([1:3,12:end])],'roi')
    end
end


end
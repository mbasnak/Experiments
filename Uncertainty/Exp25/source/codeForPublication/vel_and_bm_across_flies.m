%Code to analyze the relationship between the fly's velocity and the bump
%magnitude across flies

clear all; close all;

%% Load data

path = uigetdir('Z:\Wilson Lab\Mel\Experiments\Uncertainty\Exp25\data\Experimental');
%Get folder names
folderContents = dir(path);

%Load the summary data of the folder that correspond to experimental flies
for content = 1:length(folderContents)
   if contains(folderContents(content).name,'60D05')
       data(content) = load([folderContents(content).folder,'\',folderContents(content).name,'\analysis\vel_bm_data.mat']);
   end
end
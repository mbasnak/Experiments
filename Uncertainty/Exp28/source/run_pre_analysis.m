%Function to run the ROI_analysis and the PB analysis on all imaging
%sessions.

function run_pre_analysis(parentDir)

%Get the folder's contents
folder_contents = dir([parentDir,'\2p']);

%Retain the session ids
sessions = [];
for content = 1:length(folder_contents)
    if contains(folder_contents(content).name,'sid')
        if contains(folder_contents(content).name(1:3), 'sid')
            sessions(content) = str2num(strtok(folder_contents(content).name(5:6),'_'));
        else
            sessions(content) = NaN;
        end
    else
        sessions(content) = NaN;
    end
end
%Clean        
sessions = sort(sessions(~isnan(sessions)));

%Run the pre-analysis for all the sessions in this fly
for session = 1:length(sessions)
   %ROI_analysis(parentDir,sessions(session),0);
   PB_data_analysis(parentDir,sessions(session),0);
end

end
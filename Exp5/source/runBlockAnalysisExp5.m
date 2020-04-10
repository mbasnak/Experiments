clear all; close all;

% Navigate to the experiment data folder
cd 'Z:\Wilson Lab\Mel\Experiments\Exp5\data'
path = 'Z:\Wilson Lab\Mel\Experiments\Exp5\data\';

% List the dates in which the experiment was run
dates = dir();


% For all the actual dates
for i = 1:length(dates)
    if strfind(dates(i).name,'201') ~=0
        
        % Move to that date's folder
        newPath = [path,dates(i).name,'\experimental flies'];
        cd(newPath)
        
        % List the flies for that day
        flies = dir();
        
        % For all the actual flies
        for j = 1:length(flies)
            if strfind(flies(j).name,'flyNum') ~=0
           
                % Move to that fly's folder
                NewPath = [newPath,'\',flies(j).name,'\dataFromAnalysis'];
                cd(NewPath)
        
                blockAnalysis(NewPath)
                
            end
        end               
          
    end
    
end
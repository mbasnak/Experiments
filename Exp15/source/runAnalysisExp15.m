% Run the analysis of experiment 15 for every trial for every fly and save the plots

clear all; close all;

% Navigate to the experiment data folder
cd 'Z:\Wilson Lab\Mel\Experiments\Exp15\data'
path = 'Z:\Wilson Lab\Mel\Experiments\Exp15\data\';

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
                NewPath = [newPath,'\',flies(j).name,'\'];
                cd(NewPath)
        
                % List the trials for that fly
                trials = dir();
            
                % For all the actual trials
                for t = 1:length(trials)
                    if (strfind(trials(t).name,'Hallway') ~=0) 
                                  
                            singleTrialAnalysisExp15(NewPath,trials(t).name);
                   
                    end
                end
            end
        end               
          
    end
    
end
% Run the analysis of experiment 1 for every trial for every fly and save the plots

clear all; close all;

% Navigate to the experiment data folder
cd 'Z:\Wilson Lab\Mel\Experiments\Exp4\data'
path = 'Z:\Wilson Lab\Mel\Experiments\Exp4\data\';

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
                    if strfind(trials(t).name,'dataE') ~=0
                
                    % Run the analysis for the bar jump data                 
                    singleTrialAnalysisExp4(NewPath,trials(t).name); 
                    
                    elseif strfind(trials(t).name,'EmptyT') ~=0
                    
                    % Run the analysis for the empty trial data                 
                    EmptyTrialAnalysisExp4(NewPath,trials(t).name); 
                    
                    end
                end
            end
        end               
          
    end
    
end
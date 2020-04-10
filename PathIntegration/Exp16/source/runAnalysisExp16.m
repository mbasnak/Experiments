% Run the analysis of experiment 16 for every trial for every fly and save the plots

clear all; close all;

% Navigate to the experiment data folder
cd 'Z:\Wilson Lab\Mel\Experiments\PathIntegration\Exp16\data\experimental flies'
path = 'Z:\Wilson Lab\Mel\Experiments\PathIntegration\Exp16\data\experimental flies';

% List the genotypes
genotypes = dir();


% For all the actual genotypes
for i = 1:length(genotypes)
    if strfind(genotypes(i).name,'o') ~=0 %since they all include the letter 'o' in their name, use that
        
        % Move to that date's folder
        newPath = [path,'\', genotypes(i).name];
        cd(newPath)
        
        % List the flies for that genotpe
        flies = dir();
        
        % For all the actual flies
        for j = 1:length(flies)
            if strfind(flies(j).name,'20') ~=0
           
                % Move to that fly's folder
                NewPath = [newPath,'\',flies(j).name,'\ball'];
                cd(NewPath)
                
                % Run the analysis      
                singleTrialAnalysisExp16(NewPath)                   

            end
        end               
          
    end
    
end







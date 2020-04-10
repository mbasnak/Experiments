% Run the analysis of experiment 16 for every trial for every fly and save the plots

clear all; close all;

% Navigate to the experiment data folder
cd 'Z:\Wilson Lab\Mel\Experiments\PathIntegration\Exp17\data\experimental flies'
path = 'Z:\Wilson Lab\Mel\Experiments\PathIntegration\Exp17\data\experimental flies';

% List the genotypes
genotypes = dir();


% For all the actual genotypes
for i = 1:length(genotypes)
    if strfind(genotypes(i).name,'r') ~=0 %since they all include the letter 'r' in their name, use that
        
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
                
                % Find the files in that folder to define the relevant
                % session
                files = dir(NewPath);
                for j = 1:length(files)
                    if strfind(files(j).name,'bdata') ~=0
                        usefulfiles{j} = files(j).name;
                        sid(j) = str2num(usefulfiles{j}(55:55));
                    end
                end
                session = mode(sid);
                
                % Run the analysis   
                singleTrialAnalysisExp17(NewPath,session)                   

            end
        end               
          
    end
    
end







clear all; close all;

%% Fill out struct

%Initialize empty struct
allDataToAnalyze = struct();

%Add the parent directory
allDataToAnalyze(1).parentDir = 'Z:\Wilson Lab\Mel\Experiments\Uncertainty\Exp30\data\high_contrast\20210304_60D05_7f_fly4';
allDataToAnalyze(2).parentDir = 'Z:\Wilson Lab\Mel\Experiments\Uncertainty\Exp30\data\high_contrast\20210305_60D05_7f_fly4';
allDataToAnalyze(3).parentDir = 'Z:\Wilson Lab\Mel\Experiments\Uncertainty\Exp30\data\high_contrast\20210309_60D05_7f_fly2';
allDataToAnalyze(4).parentDir = 'Z:\Wilson Lab\Mel\Experiments\Uncertainty\Exp30\data\high_contrast\20210311_60D05_7f';
allDataToAnalyze(5).parentDir = 'Z:\Wilson Lab\Mel\Experiments\Uncertainty\Exp30\data\high_contrast\20210311_60D05_7f_fly3';
allDataToAnalyze(6).parentDir = 'Z:\Wilson Lab\Mel\Experiments\Uncertainty\Exp30\data\high_contrast\20210312_60D05_7f';
allDataToAnalyze(7).parentDir = 'Z:\Wilson Lab\Mel\Experiments\Uncertainty\Exp30\data\high_contrast\20210312_60D05_7f_fly2';
allDataToAnalyze(8).parentDir = 'Z:\Wilson Lab\Mel\Experiments\Uncertainty\Exp30\data\high_contrast\20210312_60D05_7f_fly3';


%Add the session ids
allDataToAnalyze(1).sids = 1;
allDataToAnalyze(2).sids = [0,1];
allDataToAnalyze(3).sids = [2,3];
allDataToAnalyze(4).sids = [2,3];
allDataToAnalyze(5).sids = [4,8];
allDataToAnalyze(6).sids = [1,2,3];
allDataToAnalyze(7).sids = [1,2];
allDataToAnalyze(8).sids = [0,1];


%Add the trial ids
allDataToAnalyze(1).tids = 0;
allDataToAnalyze(2).tids = [0,0];
allDataToAnalyze(3).tids = [0,0];
allDataToAnalyze(4).tids = [0,0];
allDataToAnalyze(5).tids = [0,0];
allDataToAnalyze(6).tids = [0,0,0];
allDataToAnalyze(7).tids = [0,0];
allDataToAnalyze(8).tids = [0,0];

% 
% %% Run PB analysis for all of these
% 
% for fly = 1:length(allDataToAnalyze)
%     for session = 1:length(allDataToAnalyze(fly).sids)
%         PB_data_analysis(allDataToAnalyze(fly).parentDir,allDataToAnalyze(fly).sids(session),0)
%     end
% end

%% Run bar jump analysis for all of these

for fly = 1:length(allDataToAnalyze)
    for session = length(allDataToAnalyze(fly).sids) %the bar jump session is the last session
        run_bar_jump_analysis(allDataToAnalyze(fly).parentDir,allDataToAnalyze(fly).sids(session))
    end
end
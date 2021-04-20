%Analysis for experiment #25: 3 block experiment where the first block is
%darkness/blue bar in closed loop/darkness, the second block is a bar
%changing contrasts and jumping and the third block is a horizontal bar
%changing contrasts

clear all; close all;

%get the path
path = uigetdir('Z:\Wilson Lab\Mel\Experiments\Uncertainty\Exp25\data');

%List data
analysisFiles = dir([path,'\analysis']);
%remove rows that don't actually have files
for row = 1:length(analysisFiles)
    if ~contains(analysisFiles(row).name,'analysis')
        analysisFiles(row) = [];
    end
end

%% Analyze the first block

%Load data for this block





%% Analyze the second block

clear all; close all

%Load data
[fileName,path] = uigetfile();
load([path,fileName]);





%% Analyze the third block

%identify changes in stim (ypanels)

%look at polar histograms of heading for different intensities

%plot median circ heading vs bar contrast

%plot heatmap of bump across experiment

%compute bump mag and with at half max throughout

%compare bump width at half max and bump mag vs bar contrast

%polar histogram for offset in different conditions

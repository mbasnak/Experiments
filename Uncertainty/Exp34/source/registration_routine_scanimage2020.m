%this loads a single trial data, registers it, and outputs the processed
%and refernce images for the ROI selection - made by JL

function registration_routine_scanimage2020(parentDir, sid, tid)
%===================================================================================================
%% Load a single trial data, initial process, and registration.
%% Modified by JL 1/17/2018, commented by MB 20190825
% 
% INPUTS:
%       parentDir   = directory containing the raw data that you want to process.
%
%       sid        = session id
%
%       tid         = trial id
%===================================================================================================

%% Declare variables
tic;

import ScanImageTiffReader.ScanImageTiffReader; %imports functinos necessary to open the scanimage tiff files

global slash;
if isunix() == 1 %if running in Linux or mac
    slash = '/'; %define the slash this way
else %if running in windows
    slash = '\'; %define the other way
end
set(0,'DefaultTextInterpreter','none');

% Identify the data file for the session
imagingDir = [parentDir slash '2p' slash];
ballDir = [parentDir slash 'ball' slash];

expression = ['*sid_' num2str(sid) '_tid_' num2str(tid) '_*'];
imagingFile = dir(fullfile(imagingDir, expression));
expression2 = ['bdata' '*sid_' num2str(sid) '_tid_' num2str(tid) '.mat'];
ballFile = dir(fullfile(ballDir, expression2));

disp(ballDir)
ballFile

expDate = parentDir(end-7:end);%I'll have to modifiy this because it's not the same in mine

% Make session folder for new files if necessary
sessionDir = [imagingDir sprintf('sid_%.0f_tid_%.0f',sid,tid)];
if ~isdir(sessionDir)
    mkdir(sessionDir);
end
%This makes a folder for each trial inside the 2p directory

%% Load the .tif file
% We want metadata: nFrames, numLines, numPixels, nPlanes, nVolumes,
% nChannels, reshaped.
toc;tic;

reader = ScanImageTiffReader(fullfile(imagingDir, imagingFile(1).name));
rawFile_original = reader.data();
metadata_char = reader.metadata();
%descriptions = reader.descriptions();
list = strsplit(metadata_char,'RoiGroups');%,'CollapseDelimiters',true);
list = char(list{1});
list = list(1:end-7);
eval(list); %I don't understand what this is doing. we're not storing any output, but I think it might be creating the struct SI
% This should produce a class called SI with all the metadata.
%[rawFile, metadata, SI_metadata] = read_tif(fullfile(imagingDir, imagingFile(1).name));

%% Reshape image
rawFile_original = permute(rawFile_original,[2 1 3]); %permute 1st and 2nd dimensions of the image data
dimensions = size(rawFile_original); %store imaging dimensions
if SI.hChannels.channelSave(end) == 2 %if we're saving two channels
    chan1 = rawFile_original(:,:,1:2:end);
    chan2 = rawFile_original(:,:,2:2:end);
    rawFile_channels(1:dimensions(1), 1:dimensions(2), 1:dimensions(3)/2) = chan1;
    rawFile_channels(1:dimensions(1), 1:dimensions(2), dimensions(3)/2+1:dimensions(3)) = chan2;
else
    rawFile_channels = rawFile_original;
end
rawFile = reshape(rawFile_channels, dimensions(1), dimensions(2), SI.hStackManager.numFramesPerVolume, SI.hStackManager.numVolumes, SI.hChannels.channelSave(end));
%we're reshaping the image file so that it now has 4 dimensions: the first
%two are the x and y coordinates of an imaging plane, the third one is the
%number of frames per volume, and the last one is the number of volumes
%taken in the trial.

%% Calculate the number of flyback frames
toc;tic;

ballData = load(fullfile(ballDir, ballFile.name));
bdata = ballData.trial_bdata;
samples = size(bdata, 1); %total number of frames acquired, taken from the behavioral data
index = floor(samples/SI.hStackManager.numVolumes); %index gives you the number of behavioral frames per imaging volume acquired.
second_sweep = bdata(index:index*2+1,end);%This saves the data for the piezo mvt for a volume acquired.
[~,max_index] = max(second_sweep);
[~,min_index] = min(second_sweep);
difference = abs(max_index - min_index);
fraction = difference/index;
if fraction > 0.5
    flyback = 4; % just in case this is totally wrong...
else
    flyback = round(fraction*SI.hStackManager.numFramesPerVolume); %how many flyback frames
end

%% Deal with pre-processing
toc;tic;

nChannels = SI.hChannels.channelSave(end);

if nChannels == 1
    chanData_1 = squeeze(rawFile(:,:,flyback+1:end,:,1)); %for each volume take out the flyback frames
    minimum = min(chanData_1(:)); 
    chanData_1 = chanData_1-minimum+1; %it seems like you do this to have every value be positive
else
    chanData_1 = squeeze(rawFile(:,:,flyback+1:end,:,1));
    minimum1 = min(chanData_1(:)); 
    chanData_1 = chanData_1-minimum1+1;
    chanData_2 = squeeze(rawFile(:,:,flyback+1:end,:,2));
    minimum2 = min(chanData_2(:)); 
    chanData_2 = chanData_2-minimum2+1;
end

%% Save reference image
toc;tic;

disp('Saving reference mean images...')
if nChannels > 1
    % Use the red channel to create and save reference images
    channelNum = 2;
    refImages = [];
    refImages = squeeze(mean(chanData_2, 4));
    save(fullfile(sessionDir, sprintf('sid_%.0f_tid_%.0f_refImages.mat', sid, tid)), 'refImages', 'channelNum', 'flyback');
else
    % Use the GCaMP channel to create and save reference images
    channelNum = 1;
    refImages = [];
    refImages = squeeze(mean(chanData_1, 4));
    save(fullfile(sessionDir, sprintf('sid_%.0f_tid_%.0f_refImages.mat', sid, tid)), 'refImages', 'channelNum', 'flyback');
end

%% Collect data for SAVING.
toc;tic;
% Save scanimage data from the first trial
scanImageInfo = SI;

% Get file name
fName = imagingFile.name;

%% Rigid motion registration, and SAVE DATA
% to save: regProduct, trType, fName, expDate, scanImageInfo
toc;tic;

if nChannels == 1
    % registration
    regProduct = normcorre_registration(chanData_1); %run the registration
    % save file
    channelNum = 1;
    fileName = ['sid_', num2str(sid), '_tid_', num2str(tid), '_Chan_1_sessionFile.mat'];
    savefast(fullfile(sessionDir, ['rigid_', fileName]), 'regProduct', 'fName', 'expDate', 'scanImageInfo', 'flyback');
    refImages = squeeze(mean(regProduct, 4));
    save(fullfile(sessionDir, sprintf('sid_%.0f_tid_%.0f_refImages_reg_channel_%.0f.mat', sid, tid, channelNum)), 'refImages', 'channelNum', 'flyback');
else
    % channel 1
    channelNum = 1;
    regProduct = normcorre_registration(chanData_1);
    fileName1 = ['sid_', num2str(sid), '_tid_', num2str(tid), '_Chan_1_sessionFile.mat'];
    savefast(fullfile(sessionDir, ['rigid_', fileName1]), 'regProduct', 'fName', 'expDate', 'scanImageInfo', 'flyback');
    refImages = squeeze(mean(regProduct, 4));
    save(fullfile(sessionDir, sprintf('sid_%.0f_tid_%.0f_refImages_reg_channel_%.0f.mat', sid, tid, channelNum)), 'refImages', 'channelNum', 'flyback');
    
    % channel 2
    channelNum = 2;
    regProduct = normcorre_registration(chanData_2);
    fileName2 = ['sid_', num2str(sid), '_tid_', num2str(tid), '_Chan_2_sessionFile.mat'];
    savefast(fullfile(sessionDir, ['rigid_', fileName2]), 'regProduct', 'fName', 'expDate', 'scanImageInfo', 'flyback');
    refImages = squeeze(mean(regProduct, 4));
    save(fullfile(sessionDir, sprintf('sid_%.0f_tid_%.0f_refImages_reg_channel_%.0f.mat', sid, tid, channelNum)), 'refImages', 'channelNum', 'flyback');
end

toc;
end
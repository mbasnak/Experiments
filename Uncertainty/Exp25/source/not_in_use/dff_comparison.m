%Code to compare the dff obtained with different methods


%% Load data

%clear all; close all;

[file,path] = uigetfile('Z:\Wilson Lab\Mel\Experiments\Uncertainty\Exp25\data\Experimental');
load([path,file])

%determine the sid we're working with, to save the plots to specific folder
%later
sid = file(10:14);


%% Identify changes in stim (darkness to stim)

%Plot the signal from the y dimension of the panels
figure,
subplot 121,
plot(data.fr_y_ds)
title('Panels y signal');
xlabel('Time (downsampled frames)');
ylabel('Voltage');
xlim([0 length(data.fr_y_ds)]);
subplot 122,
plot(abs(diff(data.fr_y_ds)))
title('Frames where the stimulus changes');
ylim([-1 6]);
xlim([0 length(data.fr_y_ds)]);
xlabel('Time (downsampled frames)');


%Obtain the frames where the stim changes using the derivative of the
%signal from the y dimension of the panels
changeStim = find(abs(diff(data.fr_y_ds))>0.2);
changeContrast = changeStim(2:2:end); %every other jump is a change in contrast

%% Identify jumps (using the ypanels signal)

%The jumps are the odd changes in y dim voltage
jumps = changeStim(1:2:end);

%Recover the y dimensions using the voltage
ydimensions = 24;
maxYVoltage = 10;
yDim = data.fr_y_ds*ydimensions/maxYVoltage;

figure,
plot(yDim)
xlim([0 length(yDim)]);
title('y dimension');
ylabel('Y dim'); xlabel('Time');

%Use the pos function number to define the direction of the jumps
%load the run object
run_obj_dir = dir([path(1:end-9),'ball\runobj']);
for fileName = 1:length(run_obj_dir)
    if (contains(run_obj_dir(fileName).name,file(10:14)) == 1)
        run_obj_file_name = run_obj_dir(fileName).name;
        run_obj_folder_name = run_obj_dir(fileName).folder;
    end
end
load([run_obj_folder_name,'\',run_obj_file_name])
pos_function = run_obj.function_number;

%Assign the bar jumps and the contrast changes based on the pos function
%The vector 'ContrastOrder' reflects the order in which the intensities
%appear such that ContrastOrder(1) shows in which position of the
%experiment the bar of intensity 1 appears.
%The vector 'Intensities' reflects the intensity that appears in each
%position, such that Intensities(1) returns the intensity level of the bar
%that went first in the experiment
if pos_function == 190   
    Jumps = [90, -90, 90, 90, -90, -90];
    ContrastOrder = [2,6,5,3,1,4]; 
    Intensities = [5,1,4,6,3,2];
elseif pos_function == 185  %This is for the fly from 8/28 with the previous protocol
    Jumps = [-90, 90, -90, 90, 90, 180];
    ContrastOrder = [6,4,3,1,5,2];
    Intensities = [4,6,3,2,5,1];
elseif pos_function == 191   
    Jumps = [90, -90, -90, -90, 90, 90]; 
    ContrastOrder = [5,3,6,2,1,4];
    Intensities = [5,4,2,6,1,3];
elseif pos_function == 192 
    Jumps = [-90, -90, -90, 90, 90, 90];
    ContrastOrder = [1,4,6,3,2,5];
    Intensities = [1,5,4,2,6,3];
else
    Jumps = [90, -90, 90, 90, -90, -90];
    ContrastOrder = [2,5,1,4,3,6];
end


%% Heatmap with bar jumps and stim changes

% Plot original
figure,
subplot(3,2,1)
imagesc(original_data.dff_matrix)
colormap(gray)
xlabel('Time (frames)');
ylabel('PB glomerulus');
title('Using the fifth percentile as baseline');

%using the 10th percentile
subplot(3,2,3)
imagesc(original_data.dff_matrix_10)
colormap(gray)
xlabel('Time (frames)');
ylabel('PB glomerulus');
title('Using the tenth percentile as baseline');

%without bad frames
subplot(3,2,5)
imagesc(original_data.dff_matrix_all_frames)
colormap(gray)
xlabel('Time (frames)');
ylabel('PB glomerulus');
title('Using the tenth percentile as baseline');

% Plot selecting only certain slices
subplot(3,2,2)
imagesc(data.dff_matrix)
colormap(gray)
xlabel('Time (frames)');
ylabel('PB glomerulus');
title('Using the fifth percentile as baseline');

%using the 10th percentile
subplot(3,2,4)
imagesc(data.dff_matrix_10)
colormap(gray)
xlabel('Time (frames)');
ylabel('PB glomerulus');
title('Using the tenth percentile as baseline');

%without bad frames
subplot(3,2,6)
imagesc(data.dff_matrix_all_frames)
colormap(gray)
xlabel('Time (frames)');
ylabel('PB glomerulus');
title('Using the tenth percentile as baseline');


%% Difference between original and other

figure,
subplot(5,1,1)
imagesc(original_data.dff_matrix-original_data.dff_matrix_10)
colorbar

subplot(5,1,2)
imagesc(original_data.dff_matrix-original_data.dff_matrix_all_frames)
colorbar

subplot(5,1,3)
imagesc(original_data.dff_matrix-data.dff_matrix)
colorbar

subplot(5,1,4)
imagesc(original_data.dff_matrix-data.dff_matrix_10)
colorbar

subplot(5,1,5)
imagesc(original_data.dff_matrix-data.dff_matrix_all_frames)
colorbar


%There doesn't seem to be much difference

%% Looking at the offset

figure,

subplot(6,1,1)
phase = rad2deg(original_data.phase);
offset = wrapTo180(-correctedHeading-phase);
changeOffset = abs([0,diff(offset)]);
offsetToPlot = offset;
offsetToPlot(changeOffset>100==1) = NaN;
plot(data.time,offsetToPlot,'LineWidth',1.5,'color','k')
%add the changes in stim
jumpTimes = data.time(jumps);
for change = 1:length(changeContrast)
   line([data.time(changeContrast(change)) data.time(changeContrast(change))], [-180 180], 'LineWidth', 2, 'color', 'r'); 
end
%add the jumps
for jump = 1:length(jumps)
   line([jumpTimes(jump) jumpTimes(jump)], [-180 180], 'color', 'b','LineStyle','--', 'LineWidth', 2); 
end
ylim([-180 180]);
ylabel('Deg'); xlabel('Time (sec)');
title('Standard offset');

subplot(6,1,2)
phase = rad2deg(original_data.phase_10);
offset = wrapTo180(-correctedHeading-phase);
changeOffset = abs([0,diff(offset)]);
offsetToPlot = offset;
offsetToPlot(changeOffset>100==1) = NaN;
plot(data.time,offsetToPlot,'LineWidth',1.5,'color','k')
%add the changes in stim
jumpTimes = data.time(jumps);
for change = 1:length(changeContrast)
   line([data.time(changeContrast(change)) data.time(changeContrast(change))], [-180 180], 'LineWidth', 2, 'color', 'r'); 
end
%add the jumps
for jump = 1:length(jumps)
   line([jumpTimes(jump) jumpTimes(jump)], [-180 180], 'color', 'b','LineStyle','--', 'LineWidth', 2); 
end
ylim([-180 180]);
ylabel('Deg'); xlabel('Time (sec)');
title('Offset with baseline 10th percentile');

subplot(6,1,3)
phase = rad2deg(original_data.phase_all_frames);
offset = wrapTo180(-correctedHeading-phase);
changeOffset = abs([0,diff(offset)]);
offsetToPlot = offset;
offsetToPlot(changeOffset>100==1) = NaN;
plot(data.time,offsetToPlot,'LineWidth',1.5,'color','k')
%add the changes in stim
jumpTimes = data.time(jumps);
for change = 1:length(changeContrast)
   line([data.time(changeContrast(change)) data.time(changeContrast(change))], [-180 180], 'LineWidth', 2, 'color', 'r'); 
end
%add the jumps
for jump = 1:length(jumps)
   line([jumpTimes(jump) jumpTimes(jump)], [-180 180], 'color', 'b','LineStyle','--', 'LineWidth', 2); 
end
ylim([-180 180]);
ylabel('Deg'); xlabel('Time (sec)');
title('Offset without removing bad frames');

subplot(6,1,4)
phase = rad2deg(data.phase);
offset = wrapTo180(-correctedHeading-phase);
changeOffset = abs([0,diff(offset)]);
offsetToPlot = offset;
offsetToPlot(changeOffset>100==1) = NaN;
plot(data.time,offsetToPlot,'LineWidth',1.5,'color','k')
%add the changes in stim
jumpTimes = data.time(jumps);
for change = 1:length(changeContrast)
   line([data.time(changeContrast(change)) data.time(changeContrast(change))], [-180 180], 'LineWidth', 2, 'color', 'r'); 
end
%add the jumps
for jump = 1:length(jumps)
   line([jumpTimes(jump) jumpTimes(jump)], [-180 180], 'color', 'b','LineStyle','--', 'LineWidth', 2); 
end
ylim([-180 180]);
ylabel('Deg'); xlabel('Time (sec)');
title('Standard offset with select z slices');

subplot(6,1,5)
phase = rad2deg(data.phase_10);
offset = wrapTo180(-correctedHeading-phase);
changeOffset = abs([0,diff(offset)]);
offsetToPlot = offset;
offsetToPlot(changeOffset>100==1) = NaN;
plot(data.time,offsetToPlot,'LineWidth',1.5,'color','k')
%add the changes in stim
jumpTimes = data.time(jumps);
for change = 1:length(changeContrast)
   line([data.time(changeContrast(change)) data.time(changeContrast(change))], [-180 180], 'LineWidth', 2, 'color', 'r'); 
end
%add the jumps
for jump = 1:length(jumps)
   line([jumpTimes(jump) jumpTimes(jump)], [-180 180], 'color', 'b','LineStyle','--', 'LineWidth', 2); 
end
ylim([-180 180]);
ylabel('Deg'); xlabel('Time (sec)');
title('Offset with select z slices and 10th percentile baselines');

subplot(6,1,6)
phase = rad2deg(data.phase_all_frames);
offset = wrapTo180(-correctedHeading-phase);
changeOffset = abs([0,diff(offset)]);
offsetToPlot = offset;
offsetToPlot(changeOffset>100==1) = NaN;
plot(data.time,offsetToPlot,'LineWidth',1.5,'color','k')
%add the changes in stim
jumpTimes = data.time(jumps);
for change = 1:length(changeContrast)
   line([data.time(changeContrast(change)) data.time(changeContrast(change))], [-180 180], 'LineWidth', 2, 'color', 'r'); 
end
%add the jumps
for jump = 1:length(jumps)
   line([jumpTimes(jump) jumpTimes(jump)], [-180 180], 'color', 'b','LineStyle','--', 'LineWidth', 2); 
end
ylim([-180 180]);
ylabel('Deg'); xlabel('Time (sec)');
title('Offset with select z slices and without removing the bad frames');
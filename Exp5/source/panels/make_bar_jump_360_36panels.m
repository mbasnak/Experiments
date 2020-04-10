% Make function to make a bar jump

%Create the function
PositionArray = randperm(96); %create a vector with random positions from 1 to 96
PatternSpeed = 5; %This is the speed of the pattern in Hz
NumSecStillPattern = 200; %How many seconds we want the bar not to jump for
NumFramesStillPattern = PatternSpeed*NumSecStillPattern; %number of frames where the pattern won't jump (5000 to make it 100 s at 50 Hz)
Position = repmat(PositionArray,NumFramesStillPattern,1); %Create the position function by repeating each position for 100 s
func = reshape(Position,[1,(size(Position,1))*(size(Position,2))]); %Reshape to get the vector we need

%Visualize it
time = linspace(0,NumSecStillPattern*96,size(func,2));
figure, plot(time,func)
xlabel('Time (sec)'); ylabel('Y dimension read from the pattern');
title('Function to make the bar jump');

save('200secfunction5Hz.mat','time','func');

%% SAVE position function place to be put on the SD card:
% place to save patterns to be put on the SD card:
 directory_name = 'C:\Users\Wilson\Desktop\panels-matlab_071618\functions\mel360_36panels';
 str_x = [directory_name '\position_function_bar_jump_every_' num2str(NumSecStillPattern) 'sec_5Hz']; 
 save(str_x, 'func'); 
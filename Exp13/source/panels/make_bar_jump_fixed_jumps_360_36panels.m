% Make function to make a bar jump fixed degrees

%Create the function
%Positions = [50,26,2,73,2,26,50,73,2,73,2,26,50,73,50,26,2,26,50,73,2,26,2,73,2];
%PositionArray = repmat(Positions,1,2);

%PositionArray = [2,73,2,73,50,73,2,73,2,26,2,73,2,26,2,73,50,26,50,26,50,73,2,26,2,73,2,73,50,73,2,73,2,26,2,73,2,26,2,73,50,26,50,26,50,73,2,26,2,73,2,73,50,73,2,73,2,26,2,73,2,26,2,73,50,26,50,26,50,73,2,26,2,73,2,73,50,73,2,73,2,26,2,73,2,26,2,73,50,26,50,26,50,73,2,26,2];
PositionArray = [2,73,2,73,50,73,2,73,2];


PatternSpeed = 5; %This is the speed of the pattern in Hz
NumSecStillPattern = 120; %How many seconds we want the bar not to jump for
NumFramesStillPattern = PatternSpeed*NumSecStillPattern; %number of frames where the pattern won't jump (5000 to make it 100 s at 50 Hz)
Position = repmat(PositionArray,NumFramesStillPattern,1); %Create the position function by repeating each position for 100 s
func = reshape(Position,[1,(size(Position,1))*(size(Position,2))]); %Reshape to get the vector we need

%This can be commented out, but you can use it to add a middle block in
%which y reaches a second set of dimensions (97:192).
%func(4501:4501+4499)=func(1:4500)+95;
%func(9001:9001+4499)=func(1:4500);


%Visualize it
time = linspace(0,NumSecStillPattern*96,size(func,2));
figure, plot(func)
xlabel('frames'); ylabel('Y dimension read from the pattern');
title('Function to make the bar jump');


%% SAVE position function place to be put on the SD card:
% place to save patterns to be put on the SD card:
 directory_name = 'C:\Users\Wilson\Desktop\panels-matlab_071618\functions\mel360_36panels';
 str_x = [directory_name '\position_function_166_bar_jump90deg_every_' num2str(NumSecStillPattern) 'sec_5Hz5']; 
 save(str_x, 'func'); 
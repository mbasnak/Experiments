%%  Moving stimulus to one side
% This will move the stimulus to a single side

clear all;

%% Parameters to set for each stimulus

PANELS_FRAME_RATE = 50; %Hz
%chr = 'left';
chr = 'right';

%% 

POSITION_FUNCTION_LENGTH = 1000; % this how many frames long these normally are... set by panels
% does this not depend on how long you want your stim to be for?
% Or how many LEDs you have?
numOfPanelsAcross = 12;
numOfPanelsVertically = 3;
LEDdotsPerPanel = 8; 
LEDdotsAcross = numOfPanelsAcross * LEDdotsPerPanel; 

PATTERN_SPEED_DEG_PER_SEC = 60;% deg/s  this seems like a reasonable speed for the stimulus to move at
DEGREES_PRE_PIXEL = 360/LEDdotsAcross; % ~deg

XDimMin = 0;% 
XDimMax = LEDdotsAcross;% 

PositionArray = []; % initialize an empty position array

pixelPerSecond = PATTERN_SPEED_DEG_PER_SEC / DEGREES_PRE_PIXEL; % corresponding LED pixel per second speed

framesDwellPerPixel =  round(PANELS_FRAME_RATE / pixelPerSecond); % frames/pixel=( (frames/s) / (pixel/s) )
% how many frames should be spent at each function position.
% you need to round it for it to work

% Determine whether the stimulus is moving right or left
if isequal(chr, 'left')

currFrameCounter = XDimMin; % start the counter at 0
while (length (PositionArray) < POSITION_FUNCTION_LENGTH) && (currFrameCounter <= XDimMax)
    
    addToArray = currFrameCounter * ones( 1,  framesDwellPerPixel );

    PositionArray = [PositionArray addToArray ];
    
    currFrameCounter = currFrameCounter + 1;   % move one to the left  
end
else
currFrameCounter = XDimMax - 1;
while (length (PositionArray) < POSITION_FUNCTION_LENGTH) && (currFrameCounter > XDimMin)
    
    addToArray = currFrameCounter * ones( 1,  framesDwellPerPixel );

    PositionArray = [PositionArray addToArray ];
    
    currFrameCounter = currFrameCounter - 1;   % move one to the right
end

end


NUM_TIMES_TO_REPEAT_GRATING_ROTAION = 2;
func = [ repmat(PositionArray , 1, NUM_TIMES_TO_REPEAT_GRATING_ROTAION)];

%% SAVE position function place to be put on the SD card:
% place to save patterns to be put on the SD card:
 directory_name = 'C:\Users\Wilson\Desktop\panels-matlab_071618\functions\mel360_36panels';
 str_x = [directory_name '\position_function_055_moving' chr '_' num2str(PATTERN_SPEED_DEG_PER_SEC) 'degS']; 
 save(str_x, 'func'); % variable must be named 'func'
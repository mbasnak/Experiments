%% starfield, with y dimensions to store positions to use for optic flow

clear all;
numOfPanelsAcross = 12;
numOfPanelsVertically = 3;
LEDdotsPerPanel = 8;  
LEDdotsAcross = numOfPanelsAcross * LEDdotsPerPanel; 
LEDdotsVertically = numOfPanelsVertically * LEDdotsPerPanel;

%Save general infomation about pattern layout
pattern.x_num = LEDdotsAcross;
pattern.y_num = 96; 	

pattern.num_panels = 36; 
pattern.gs_val = 1; 	

Pats = zeros(LEDdotsVertically, LEDdotsAcross, pattern.x_num, pattern.y_num); 	%initializes the array with zeros

% Construct the grating pattern in Y = 2, and all x values
% zeros 0 = dark, ones 1 = light
PATTERN_INDEX = 1;
gratingWidth = 3;
grating = zeros( 1 , LEDdotsAcross);

for i = 1: LEDdotsAcross  % for every one of the 96 LEDs in the x axis
    if( mod( i , gratingWidth * 2) < ( gratingWidth ) )
        % mod gives the remainder after division of i by gratingWidth*2
        grating(i) = 1; % add in light stripes at each of these locations
        % this will basically make every consecutive 3 (gratingWidth) columns light, and
        % the next 3 dark, etc...
    end
end

stripe_pattern = repmat ( grating, LEDdotsVertically,  1); % extend matrix to include the whole vertical extent of the panels

%make 2 dark ends 180 deg apart
stripe_pattern(:,[15:29,63:77])=0;
%% Make the x dimensions move circularly (this will be controlled by the animal's yaw')

% For this stimulus, the rotation occurs across every x dimension if the
% animal is yawing.

Pats(:,:,1,1) = stripe_pattern;
% 
% for x = 2:pattern.x_num % for each frame in the x axis
%     Pats(:,:,x,:) = circshift(Pats(:,:,x-1,:),[0,1]); % make the value of the frame equal to moving the 3rd dimension by one pixel
% end
% 
% %save as video to visualize moving pattern
%define colormap for visualization
map = [0 0 0.3
    0 0 0.4
    0 0 0.5
    0 0 0.6
    0 0 0.8
    0 0 1.0];

colormap(map)
% for i = 1:96
%     figure(1)  
%     imagesc(Pats(:,:,i,1))
%     title(strcat('frame',num2str(i)));
%     F(i) = getframe(gcf) ;
%     drawnow
% end
% % create the video writer with 1 fps
% writerObj = VideoWriter('Z:\Wilson Lab\Mel\FlyOnTheBall\data\Experiment15\stimVideos\rotationalFlowgratingDarkEnds.avi');
% writerObj.FrameRate = 1; % set the seconds per image
% % open the video writer
% open(writerObj);
% % write the frames to the video
% for i=1:length(F)
%     frame = F(i);        % convert the image to a frame
%     writeVideo(writerObj, frame);
% end
% % close the writer object
% close(writerObj);


%% Make the y dimension move F to B for the desired windows 

%There are 2 52 deg windows in which the animal gets front to back motion,
%although he might get a weird combination of front to back and rotational

for x = 1:pattern.x_num
    Pats(:,:,x,1) = stripe_pattern;
end
        
 for y = 2:pattern.y_num
        
        %Make the grating flow from front to back for the 2 x dimensions
        %180 deg apart from one another
        
        Pats(:,29:63,[1:14,49:61],y) = circshift(Pats(:,29:63,[1:14,49:61],y-1),[0,-1]);

        for i = 78:96
            Pats(:,i,[1:14,49:61],y) = Pats(:,i-1,[1:14,49:61],y-1);
        end
        
        for i = 2:14
            Pats(:,i,[1:14,49:61],y) = Pats(:,i-1,[1:14,49:61],y-1);
        end
        
        Pats(:,77,[1:14,49:61],y) = Pats(:,82,[1:14,49:61],y-1);
        Pats(:,1,[1:14,49:61],y) = Pats(:,96,[1:14,49:61],y-1);
        
                
        %Leave the grating static for the other x dimensions when moving in
        %y
        Pats(:,:,15:48,y) = Pats(:,:,15:48,1);  
        Pats(:,:,62:96,y) = Pats(:,:,62:96,1);
                
end 


%save as video to visualize moving pattern
colormap(map)

for i = 1:96
    figure(1)  
    imagesc(Pats(:,:,3,i))
    title(strcat('frame',num2str(i)));
    F(i) = getframe(gcf) ;
    drawnow
end
% create the video writer with 1 fps
writerObj = VideoWriter('Z:\Wilson Lab\Mel\FlyOnTheBall\data\Experiment15\stimVideos\FtoBgratingDarkEnds.avi');
writerObj.FrameRate = 1; % set the seconds per image
% open the video writer
open(writerObj);
% write the frames to the video
for i=1:length(F)
    frame = F(i);        % convert the image to a frame
    writeVideo(writerObj, frame);
end
% close the writer object
close(writerObj);


%% Create pattern map

pattern.Pats = Pats; % put data in structure 

% Make panel map for 360 arena
pattern.Panel_map = [12 8 4 11 7 3 10 6 2 9 5 1;...
                     24 20 16 23 19 15 22 18 14 21 17 13;...
                     36 32 28 35 31 27 34 30 26 33 29 25];

pattern.BitMapIndex = process_panel_map(pattern);
pattern.data = Make_pattern_vector(pattern);


 
%% SAVE pattern place to save patterns to be put on the SD card:
% place to save patterns to be put on the SD card:
directory_name = 'C:\Users\Wilson\Desktop\panels-matlab_071618\Patterns\mel360_36panels';
 str = [directory_name '\Pattern033_gratings_360_XY_opticflow180deg_50degwindow_darkends']; 	% name must begin with ‘Pattern_’
 save(str, 'pattern');
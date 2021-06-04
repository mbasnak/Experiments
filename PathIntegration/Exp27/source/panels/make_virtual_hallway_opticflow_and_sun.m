%Make pattern that has low contrast gratings with front to back
%optic flow to the sides as the animal moves forward and a constant 'Sun'
%cue to the side

clear all; close all;

pattern.x_num = 96; % number of pixels in the x axis (there are 8 pixels per panel, and 12 panels)
pattern.y_num = 92; 
pattern.num_panels = 36; % there are 36 panels
pattern.gs_val = 4; % we are setting the grayscale value at 4 for defining the intensity of the blue LEDs.
%The blue bars of the gratings will be gs_val = 1, and the sun will be
%gs_val = 15;
Pats = zeros(24, 96, pattern.x_num, pattern.y_num);

%%
%-----------------------------------------------------
% Stripe optic flow
%Because of how my LED panels are arranged, the right side of the stimulus
%goes from px 71 in the front to px 22 in the back, and the left side of
%the stimulus goes from px 70 in the front to px 23 in the back
%My stim pattern will then have a structure where I separate in sections
%what happens with pxs 1-22, 23-70 and 71-96.

% stripe_pattern{1} = [zeros(24,1),repmat([ones(24,2),zeros(24,2)],1,4),ones(24,2),zeros(24,1)                               ,zeros(24,1),repmat([ones(24,2),zeros(24,2)],1,11),ones(24,2),zeros(24,1)          ,zeros(24,1),repmat([ones(24,2),zeros(24,2)],1,6),ones(24,2),zeros(24,1)];
% stripe_pattern{2} = [zeros(24,2),repmat([ones(24,2),zeros(24,2)],1,4),ones(24,2)                                           ,repmat([ones(24,2),zeros(24,2)],1,12)                                             ,zeros(24,2),repmat([ones(24,2),zeros(24,2)],1,6),ones(24,2)];
% stripe_pattern{3} = [ones(24,1),zeros(24,2),repmat([ones(24,2),zeros(24,2)],1,4),ones(24,1)                                ,ones(24,1),zeros(24,2),repmat([ones(24,2),zeros(24,2)],1,11),ones(24,1)           ,ones(24,1),zeros(24,2),repmat([ones(24,2),zeros(24,2)],1,6),ones(24,1)];
% stripe_pattern{4} = [repmat([ones(24,2),zeros(24,2)],1,5)                                                                  ,zeros(24,2),repmat([ones(24,2),zeros(24,2)],1,11),ones(24,2)                      ,repmat([ones(24,2),zeros(24,2)],1,7)];                                                                
% 
% 
% for i = 5:48
%     stripe_pattern{i} = stripe_pattern{i-4};
% end
% 
% 
% for i = 49:70
%     stripe_pattern{i} = stripe_pattern{24};
% end
% 
% for i = 71:92
%     stripe_pattern{i} = stripe_pattern{1};
% end
% 
% for y = 1:pattern.y_num %for every y dimension
%     
%     Pats(:,:,1,y) = stripe_pattern{1,y}; %the x dim = 1 will be the stripe pattern for that dimension
%     Pats(:,53:54,1,y) = 15; %add the Sun stimulus
%     
%     for x = 2:pattern.x_num
%         Pats(:,:,x,y) = circshift(Pats(:,:,x-1,y),[0,1]); %the other x dims will be shifting one px in yaw per 1 dim shift
%     end
% 
% end
% ----------------------------------------------------------------------

% Starfield Optic Flow
% Right side front to back from px 62 to 31, Left side front to back from
% px 79 to 96 and then px 1 to 14. From 63 to 78 and from 15 to 30 are set
% to darkness.

dot_density=0.20;
% Draw the initial image at xpos=l, ypos=1

star_pattern{1}=zeros(24,32); % The right side from px 31 to 62 
% Shuffle a list of linear indices
shuffled_pixels=randperm(numel(star_pattern{1})); 
numDots=round( length(shuffled_pixels)*dot_density );
for d=1:numDots
    star_pattern{1}(shuffled_pixels(d))=1;
end

for y=2:pattern.y_num
    star_pattern{y}=circshift(star_pattern{y-1},[0,-1]);
end

for y = 1:pattern.y_num % for every y dimension
    
    % px 1-14 = px 44-31
    Pats(:,:,1,y) = [flip(star_pattern{y}(:,1:14),2),zeros(24,16),star_pattern{y},zeros(24,16),flip(star_pattern{y}(:,15:32),2)];
    Pats(:,63:64,1,y) = 15; %add the Sun stimulus
    
    for x = 2:pattern.x_num
        Pats(:,:,x,y) = circshift(Pats(:,:,x-1,y),[0,1]); %the other x dims will be shifting one px in yaw per 1 dim shift
    end

end

%%

% map = [0 0 0; 0 0 1];
% 
% video = VideoWriter('hallway.avi'); %create the video object
% open(video); %open the file for writing
% for i=1:13 %where N is the number of images
%   hAxes = gca;
%   imagesc(hAxes, Pats(:,:,1,i))
%   colormap(hAxes , map)
%   set(gca,'XTick',[], 'YTick', [])
%   saveas(gcf,['pattern',num2str(i),'.jpg'])
%   I = imread(['pattern',num2str(i),'.jpg']); %read the next image
%   writeVideo(video,I); %write the image to file
%   pause(1)
% end
% close(video); %close the file



pattern.Pats = Pats;

% Change mapping to accomodate bad '31' wz 2021/05/19
pattern.Panel_map = [36 32 28 35 25 27 34 30 26 33 29 31;...
                     24 20 16 23 13 15 22 18 14 21 17 19;...
                     12 8 4 11 1 3 10 6 2 9 5 7];

pattern.Panel_map = pattern.Panel_map;
pattern.BitMapIndex = process_panel_map(pattern);
pattern.data = Make_pattern_vector(pattern);



%% Save data

directory_name = 'C:\Users\Wilson\Desktop\panels-matlab_071618\Patterns\mel360_36panels';
str = [directory_name '\Pattern031_newMap_vhallway_with_star_opticflow1_and_sun15'];
save(str, 'pattern'); % save the file in the specified directory
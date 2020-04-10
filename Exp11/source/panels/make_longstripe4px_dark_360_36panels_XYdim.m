% To make a 4 px light bar with both x and y dimensions, to allow for bar
% jumps

pattern.x_num = 96; 
pattern.y_num = 96; %I'm making 12 dimensions in the y direction, and I can make the bar jump to 12 different positions 
pattern.num_panels = 36;
pattern.gs_val = 3; 
Pats = zeros(24, 96, pattern.x_num, pattern.y_num);

stripe_pattern = [zeros(24,4),ones(24, 92)]; 

Pats(:,:,1,2) = stripe_pattern;

for x = 2:pattern.x_num % for each frame in the x axis
    %Pats(:,:,x,:) = ShiftMatrix(Pats(:,:,x-1,:),1,'l','y'); % make the value of the frame equal to moving the 3rd dimension by one pixel
    Pats(:,:,x,:) = circshift(Pats(:,:,x-1,:),[0,1]); % make the value of the frame equal to moving the 3rd dimension by one pixel
    
    for y = 3:pattern.y_num
        %Pats(:,:,:,y) = ShiftMatrix(Pats(:,:,:,y-1),1,'l','y');
        Pats(:,:,:,y) = circshift(Pats(:,:,:,y-1),[0,1]);
        Pats(:,:,:,1) = zeros;
    end 
end

% circshift(M,[x,y]);
pattern.Pats = Pats;

pattern.panel_map = [12 8 4 11 7 3 10 6 2 9 5 1;...
                     24 20 16 23 19 15 22 18 14 21 17 13;...
                     36 32 28 35 31 27 34 30 26 33 29 25];
                 
pattern.Panel_map = pattern.panel_map;
pattern.BitMapIndex = process_panel_map(pattern);
pattern.data = Make_pattern_vector(pattern);

%% 

directory_name = 'C:\Users\Wilson\Desktop\panels-matlab_071618\Patterns\mel360_36panels';
str = [directory_name '\Pattern018_4px_darkstripe_360_36panels_XYdim'];
save(str, 'pattern');
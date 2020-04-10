% To make a 2 px light bar with both x and y dimensions, to allow for bar
% jumps

pattern.x_num = 96; 
pattern.y_num = 96; %I'm making 12 dimensions in the y direction, and I can make the bar jump to 12 different positions 
pattern.num_panels = 24;
pattern.gs_val = 1; 
Pats = zeros(16, 96, pattern.x_num, pattern.y_num);

stripe_pattern = [zeros(16, 94), ones(16, 2)]; 

Pats(:,:,1,1) = stripe_pattern;

for x = 2:pattern.x_num % for each frame in the x axis
    %Pats(:,:,x,:) = ShiftMatrix(Pats(:,:,x-1,:),1,'l','y'); % make the value of the frame equal to moving the 3rd dimension by one pixel
    Pats(:,:,x,:) = circshift(Pats(:,:,x-1,:),[0,1]); % make the value of the frame equal to moving the 3rd dimension by one pixel
    
    for y = 2:pattern.y_num
        %Pats(:,:,:,y) = ShiftMatrix(Pats(:,:,:,y-1),1,'l','y');
        Pats(:,:,:,y) = circshift(Pats(:,:,:,y-1),[0,1]);
    end 
end
% circshift(M,[x,y]);
pattern.Pats = Pats;

pattern.panel_map = [22 13 17 21 16 20 24 15 19 23 14 18; 10 1 5 9 4 8 12 3 7 11 2 6];
pattern.Panel_map = pattern.panel_map;
pattern.BitMapIndex = process_panel_map(pattern);
pattern.data = Make_pattern_vector(pattern);


directory_name = 'Z:\Wilson Lab\Mel\set up building\LED panels\Panels\panels-matlab_071618\Patterns\mel360';
str = [directory_name '\Pattern_longstripe2px_light_360_XYdim'];
save(str, 'pattern');
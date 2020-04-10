% Code to make a pattern consisting of a 2px bright blue bar on a dark
% background in a 360 deg arena
% This needs to be run in the panels-matlab folder that has the support
% functions necessary

pattern.x_num = 96; % number of pixels in the x axis (there are 8 pixels per panel, and 12 panels per row)
pattern.y_num = 1; % number of y dimensions of the stimulus; in this case, 1
pattern.num_panels = 24; % for this experiments, I am using 2 rows of LED panels in a 360 deg arena, so there are 24 panels total 
pattern.gs_val = 1; % the grayscale value was set to 1
Pats = zeros(16, 96, pattern.x_num, pattern.y_num);
% Pats file of size (L,M,N,O), where L is the number of pixel rows, 8
% per panel, M is the number of px columns, 8 per panel, N is the number of
% frames in the 'x' dimmension, and O is the number of frames in the 'y' dimension

stripe_pattern = [zeros(16, 94), ones(16, 2)]; 
% I'm initializing the light stripe pattern consists of a stripe of 16x2 pixels light over a
% 16x94 pixels dark background

Pats(:,:,1,1) = stripe_pattern;
% I'm setting the first dimension of the Pats file as the stripe_pattern

for j = 2:96 % for each frame in the x axis
    Pats(:,:,j,1) = ShiftMatrix(Pats(:,:,j-1,1),1,'l','y');
    % make the value of the frame equal to moving the 3rd dimension by one pixel
    % (i.e., circularly move the bar one pixel)
end

pattern.Pats = Pats;

pattern.panel_map = [22 13 17 21 16 20 24 15 19 23 14 18; 10 1 5 9 4 8 12 3 7 11 2 6];
% we save our panel map according to the arrangement of the LED panels in
% the board. The first 12 values are the upper LEDs

pattern.Panel_map = pattern.panel_map;
pattern.BitMapIndex = process_panel_map(pattern);
% we run the function "process_panel_map" with our pattern info.
% this returns a field in the struct pattern called BitMapIndex, that has
% for each panel number the list of x and y pixels that compose it

pattern.data = Make_pattern_vector(pattern);
% this runs the function "Make_pattern_vector" and returns the output in
% the field "data". 
% This is giving a vector of the pattern in grayscale values

directory_name = 'C:\Users\Melanie\Documents\Mel\set up building\LED panels\Panels\panels-matlab_071618\Patterns\mel360';
str = [directory_name '\Pattern_longstripe2px_light_360'];
save(str, 'pattern'); % save the file in the specified directory
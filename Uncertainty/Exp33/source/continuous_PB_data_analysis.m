function continuous_PB_data_analysis(path,sid,tid)

%Code to analyze the imaging and behavior data, with the imaging data
%analyzed using the new 'continuous' method

global slash;
if isunix() == 1 %if running in Linux or mac
    slash = '/'; %define the slash this way
else %if running in windows
    slash = '\'; %define the other way
end
set(0,'DefaultTextInterpreter','none');

%% Create the data analysis folder

data_analysis_dir = [path slash 'analysis' slash];

if(~exist(data_analysis_dir, 'dir'))
    mkdir(data_analysis_dir);
end

%% Load the imaging data

%Move to the folder of interest
cd(path)

%Load the roi data 
load(['2p' slash 'ROI' slash 'ROI_midline_sid_' num2str(sid) '_tid_' num2str(tid) '.mat']);
%Load the registered imaging stack
load(['2p' slash 'sid_' num2str(sid) '_tid_' num2str(tid) slash 'rigid_sid_' num2str(sid) '_tid_' num2str(tid) '_Chan_1_sessionFile.mat']);

%% Get the summed GCaMP7f data

summedData = squeeze(sum(regProduct,3)); %add the data across the 8 z layers making up each volume, to obtain 1 PB image per timepoint

%% Get midline coordinates and PB mask

%Find the row corresponding to the midline
for row = 1:length(roi)
    if contains(roi(row).name,'mid')
        roi_mid = row;
    end
end

%1) Pull up the values that make up the mid line along the PB
midline_coordinates = [roi(roi_mid).xi,roi(roi_mid).yi];
%2) Get the 'vector lengths' for each section of the line, to know how much
%distance of the PB they represent
midline_segment_lengths = zeros(length(midline_coordinates)-1,1);
for segment = 1:length(midline_coordinates)-1
    midline_segment_lengths(segment) = sqrt( (midline_coordinates(segment+1,1)-midline_coordinates(segment,1))^2 + (midline_coordinates(segment+1,2)-midline_coordinates(segment,2))^2 );
end
midline_distances = cumsum(midline_segment_lengths);

%3) Pull up the locations for all the px in the PB mask
PB_mask = zeros(size(regProduct,1),size(regProduct,2));
for row = 1:length(roi)
    if row~=roi_mid
        PB_mask = PB_mask + roi(row).BW;
    end
end
PB_mask(PB_mask>1) = 1;
PB_coverage = logical(PB_mask);

%% Get the normal line across each small segment of the PB midline

midV = cell(1,length(midline_coordinates)-1);
normal1 = cell(1,length(midline_coordinates)-1);
normal2 = cell(1,length(midline_coordinates)-1);
point1 = cell(1,length(midline_coordinates)-1);
point2 = cell(1,length(midline_coordinates)-1);

for segment = 1:length(midline_coordinates)-1
    
    y_part = [midline_coordinates(segment,1),midline_coordinates(segment,2)];
    x_part = [midline_coordinates(segment+1,1),midline_coordinates(segment+1,2)];
    V = y_part - x_part;
    midV{segment} = x_part + 0.5 * V;
    normal1{segment} = [ V(2), -V(1)];
    normal2{segment} = [-V(2),  V(1)];
    
    point1{segment} = [midV{segment}(1), midV{segment}(2)];
    point2{segment} = [midV{segment}(1) + normal1{segment}(1), midV{segment}(2) + normal1{segment}(2)];
end

 
%% Assign fluorescence across the PB to the corresponding closest point in the midline

%tic

midline_ff = zeros(length(midline_coordinates)-1,size(summedData,3)); %start empty midline brightness vector
PB_Image = cell(1,size(summedData,3));

for timepoint = 1:size(summedData,3)
    
    PB_Image{timepoint} = summedData(:,:,timepoint);
    PB_Image{timepoint}(PB_coverage == 0) = 0;
    
    %For each of those locations, find the distance to each normal in the
    %midline
    for row = 1:size(PB_coverage,1)
        for column = 1:size(PB_coverage,2)
            if PB_coverage(row,column) == 1
                %for every normal to the midline
                for norm_line = 1:length(normal1)
                    % Get the distance to that normal
                    dist_to_midline(norm_line) = GetPointLineDistance(column,row,point1{norm_line}(1),point1{norm_line}(2),point2{norm_line}(1),point2{norm_line}(2));
                end
                %find the minimum distance
                [~, Imin] = min(dist_to_midline);
                %add intensity value of that pixel to that point in the midline
                midline_ff(Imin,timepoint) = midline_ff(Imin,timepoint) + PB_Image{timepoint}(row,column);
            end
        end
    end
        
end

%toc

%% Get baseline and compute DF/F

baseline_f = zeros(1,length(midline_coordinates)-1);
for segment = 1:length(midline_coordinates)-1
    
    sorted_f = sort(midline_ff(segment,:));
    %Get baseline as the tenth percentile
    baseline_f(segment) = prctile(sorted_f,10);  
    
end

% Compute df/f
dff = (midline_ff'-baseline_f)./baseline_f;


%% Compute bump parameters

%tic
[bump_mag, bump_width, adj_rs, bump_pos] = fitVonMises(midline_distances,dff);
%toc

%% Load the ball data

ball_dir = [path slash 'ball' slash];
expression = ['*sid_' num2str(sid) '_tid_' num2str(tid) '*.mat'];
ball_file = dir(fullfile(ball_dir, expression));
ballData = load(fullfile(ball_dir, ball_file.name)); %load ballData
runobjFile = dir(fullfile([ball_dir,'\runobj\'], ['*_sid_' num2str(sid) '_*']));
load(fullfile([ball_dir,'\runobj\'], runobjFile.name)); %load run_obj


%% Convert the ball data panel position into an angle, bar position; fictrac position/velocities

% 1)Import behavior file
ball_file = dir(fullfile(ball_dir, expression));
ballData = load(fullfile(ball_dir, ball_file.name)); %load ballData
runobjFile = dir(fullfile([ball_dir,'\runobj\'], ['*_sid_' num2str(sid) '_*']));
load(fullfile([ball_dir,'\runobj\'], runobjFile.name)); %load run_obj
bdata_raw = ballData.trial_bdata; %get the ball data
bdata_time = ballData.trial_time; %get the trial time

% 2)Use an auxiliary function to get the different components of the behavior data
number_x = 96;
[smoothed, bdata_time_out, panel_angle, flyPosRad] = get_data_360(bdata_time, bdata_raw, number_x);

% 3)Recover relevant movement parameters
vel_for = smoothed.xVel';
vel_yaw = smoothed.angularVel;
vel_for_deg = smoothed.xVelDeg';
vel_side_deg = smoothed.yVelDeg';
total_mvt = smoothed.total_mvt;

%If you want to check what the movement data looks like, uncomment the
%following line
%plot_mvt_parameters(vel_for_deg,vel_side_deg,vel_yaw,total_mvt)

panel_y = downsample(bdata_raw(:,6), floor(4000/50));

% 4)Subsample all the variables to have the length of the number of
%volumes scanned:
volumes = length(dff); %get the volume number
time_ds = bdata_time_out(round(linspace(1, length(bdata_time_out), volumes))); 
vel_for_ds = vel_for(round(linspace(1, length(vel_for), volumes)));
vel_yaw_ds = vel_yaw(round(linspace(1, length(vel_yaw), volumes)));
vel_for_deg_ds = vel_for_deg(round(linspace(1, length(vel_for_deg), volumes)));
vel_side_deg_ds = vel_side_deg(round(linspace(1, length(vel_side_deg), volumes)));
total_mvt_ds = total_mvt(round(linspace(1, length(total_mvt), volumes)));
panel_y_ds = panel_y(round(linspace(1, length(panel_y), volumes)));

%With this convention, a positive change in panel_angle_ds implies a
%clockwise change in bar position.
panel_angle_ds = panel_angle(round(linspace(1, length(panel_angle), volumes)));

%With this convention, a positive change in flyPosRad_ds implies a
%clockwise change in heading.
flyPosRad_ds = resample(wrapToPi(flyPosRad),volumes,length(flyPosRad));


%% Offset calculation
%we will compute and define the offset as the
%circular distance between the negative of dff_pva and the heading
%(since a clockwise movement in the EB should elicit a counterclockwise movement of the fly)
offset = wrapTo360(rad2deg(circ_dist(-bump_pos',flyPosRad_ds)));


%% Save the data into the analysis folder

% General experiment info
continuous_data.parentDir = path;
continuous_data.sid = sid;
continuous_data.tid = tid;
continuous_data.time = time_ds;
continuous_data.run_obj = run_obj;
continuous_data.trial_dur = run_obj.trial_t;
continuous_data.fr_y_ds = panel_y_ds;

% Movement parameters
continuous_data.vel_for = vel_for;
continuous_data.vel_yaw = vel_yaw;
continuous_data.vel_for_ds = vel_for_ds;
continuous_data.vel_yaw_ds = vel_yaw_ds;
continuous_data.vel_side_deg_ds = vel_side_deg_ds;
continuous_data.vel_for_deg_ds = vel_for_deg_ds;
continuous_data.total_mvt_ds = total_mvt_ds;
continuous_data.heading = flyPosRad_ds;
continuous_data.heading_deg = rad2deg(flyPosRad_ds);
continuous_data.panel_angle = panel_angle_ds;

% Imaging data
continuous_data.volumes = volumes;
continuous_data.dff_matrix = dff;
continuous_data.bump_magnitude = bump_mag;
continuous_data.bump_width = bump_width;
continuous_data.adj_rs = adj_rs;
continuous_data.bump_pos = wrapToPi(bump_pos);
continuous_data.offset = offset;


% Write file
filename = [path, '\analysis\continuous_analysis_sid_' num2str(sid) '_tid_' num2str(tid) '.mat'];
save(filename, 'continuous_data');

end
function continuous_PB_data_analysis(path,sid,tid)

%Move to the folder of interest
cd(path)

%Load the roi data 
load(['2p/ROI/ROI_midline_sid_',num2str(sid),'_tid_',num2str(tid),'.mat']);
%Load the registered imaging stack
load(['2p/sid_',num2str(sid),'_tid_',num2str(tid),'/rigid_sid_',num2str(sid),'_tid_',num2str(tid),'_Chan_1_sessionFile.mat']);

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

% figure('Position',[200 50 600 800]);
% subplot(3,1,1)
% plot(midline_coordinates(:,1),midline_coordinates(:,2), 'r');
% hold on
% for segment = 1:length(midline_coordinates)-1
%     plot([midV{segment}(1), midV{segment}(1) + normal1{segment}(1)], [midV{segment}(2), midV{segment}(2) + normal1{segment}(2)], 'b');
%     plot([midV{segment}(1), midV{segment}(1) + normal2{segment}(1)], [midV{segment}(2), midV{segment}(2) + normal2{segment}(2)], 'b');
% end
% title('PB midline and normal segments');

%% Get the distance to each normal line in an example timepoint

% PB_image = summedData(:,:,1000);
% PB_image(PB_coverage == 0) = 0;
% 
% %For each of those locations, find the distance to each normal in the
% %midline
% midline_f = zeros(1,length(midline_coordinates)-1); %start empty midline brightness vector
% 
% for row = 1:size(PB_coverage,1)
%     for column = 1:size(PB_coverage,2)
%         if PB_coverage(row,column) == 1
%             %for every normal to the midline
%             dist_to_midline = zeros(length(normal1),1);
%             for norm_line = 1:length(normal1)
%                 % Get the distance to that normal
%                 dist_to_midline(norm_line) = GetPointLineDistance(column,row,point1{norm_line}(1),point1{norm_line}(2),point2{norm_line}(1),point2{norm_line}(2));
%             end
%             %find the minimum distance
%             [~ , Imin] = min(dist_to_midline);
%             %add intensity value of that pixel to that point in the midline
%             midline_f(Imin) = midline_f(Imin) + PB_image(row,column);
%         end
%     end
% end
% % 
% % Plot
% subplot(3,1,2)
% imagesc(flip(PB_image))
% title('PB activity at example timepoint');
% colormap(gray)
% 
% subplot(3,1,3)
% plot(midline_distances,midline_f)
% title('Fluorescence at example timepoint using midline');
% 
% saveas(gcf,'ExampleTimepoint.png');

% 
%% Repeat for every timepoint

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

% figure('Position',[100 100 1400 300]),
% imagesc(dff')
% colormap(flipud(gray));
% title('Dff using new method');
% saveas(gcf,'epg_activity_heatmap.png');


%% Compute bump parameters

%tic
[bump_mag, bump_width, adj_rs, u] = fitVonMises(midline_distances,dff);
%toc

%Plot
% figure,
% subplot(2,1,1)
% plot(bump_mag)
% title('Bump magnitude');
% 
% subplot(2,1,2)
% plot(bump_width)
% title('Bump width');

%% Figure out which frames to ditch

frames = [1:length(dff)];

% %Plot
% figure,
% subplot(3,1,1)
% imagesc(dff')
% colormap(flipud(gray))
% 
% subplot(3,1,2)
% plot(frames(adj_rs>=0.5),bump_mag(adj_rs>=0.5),'.r')
% hold on
% plot(frames(adj_rs<0.5),bump_mag(adj_rs<0.5),'.k')
% title('Bump magnitude');
% ylabel('DF/F');
% xlim([1 length(adj_rs)]);
% legend('adjrs >= 0.5','adjrs < 0.5');
% 
% subplot(3,1,3)
% plot(frames(adj_rs>=0.5),bump_width(adj_rs>=0.5),'.r')
% hold on
% plot(frames(adj_rs<0.5),bump_width(adj_rs<0.5),'.k')
% title('Bump witdh');
% ylabel('Radians');
% xlim([1 length(adj_rs)]);
% legend('adjrs >= 0.5','adjrs < 0.5');
% 
% %Get percentage of useful frames (good gof)
% perc_good = (sum(adj_rs>=.5)/length(adj_rs))*100;
% suptitle([num2str(round(perc_good)),'% good frames']);

%% Save data

continuous_data.dff_matrix = dff;
continuous_data.bump_magnitude = bump_mag;
continuous_data.bump_width = bump_width;
continuous_data.adj_rs = adj_rs;
continuous_data.bump_pos = u;

filename = [path, '\analysis\continuous_analysis_sid_' num2str(sid) '_tid_' num2str(tid) '.mat'];
save(filename, 'continuous_data');

end
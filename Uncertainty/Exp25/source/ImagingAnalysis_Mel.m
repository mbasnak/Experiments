%% Play the movie of activity
imgData = squeeze(sum(regProduct,3));

for frame = 1:size(imgData,3)
    
    pause_time = 0.05;
    h = imagesc(imgData(:,:,frame));
    axis equal tight
    colormap(bone)
    
    for i = 1:size(regProduct,4)
        h.CData = sum(regProduct(:,:,:,i),3);
        pause(pause_time)
    end
    
end

%% plot the maximum projection

imagesc(max(imgData,[],3)) %(max intensity over time of summed intensity over planes)
colormap(bone)
axis equal tight

%% draw a bounding box over the PB

mask = roipoly;

%% extract midline axis, and subsample for centroid locations

n_centroid = 40;
mid  = bwmorph(mask,'thin',inf);
[y,x] = find(mid);
x = smooth(x,15);
y = smooth(y,15);
idx = round(linspace(1,length(y),n_centroid));
centroids = [y(idx),x(idx)];
cmap = cbrewer2('set1',n_centroid);

clf
imagesc(max(sum(regProduct,3),[],4)) %(max intensity over time of summed intensity over planes)
colormap(bone)
hold on
plot(x,y,'r')
scatter(centroids(:,2),centroids(:,1),[],cmap,'filled')
axis equal tight

%% assign each pixel to a centroid

imgData_mask = imgData.*mask;

[y_coor,x_coor] = find(mask); %find coordinates of pixels in mask
[dists,idx] = pdist2(centroids,[y_coor,x_coor],'euclidean','smallest',1);
hold on
for i = 1:n_centroid
    scatter(x_coor(idx == i),y_coor(idx == i),'filled','MarkerFaceColor',cmap(i,:),'MarkerFaceAlpha',0.2)
end

%% Find the median activity in each group over time
imgData_2d = reshape(imgData,[],size(imgData,3));
centroid_log = false(n_centroid,size(imgData_2d,1));
for i = 1:n_centroid
    centroid_log(i, sub2ind(size(imgData),y_coor(idx==i),x_coor(idx ==i))) = true;
end

avg_intensity = centroid_log * imgData_2d ./ sum(centroid_log,2);

%% PLot both

pause_time = 0.1;

figure
subplot(2,1,1)
h(1) = plot(nan(n_centroid,1));
ylim([0,max(avg_intensity,[],'all')])

subplot(2,1,2)
h(2) = imagesc(imgData(:,:,i));
axis equal tight
colormap(bone)

for i = 1:size(avg_intensity,2)
    h(1).YData = avg_intensity(:,i);
    h(2).CData = imgData(:,:,i);
    pause(pause_time)
end
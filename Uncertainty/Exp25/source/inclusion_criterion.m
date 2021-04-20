%Code that plots the reference offset as well as the verifying offset, to
%decide whether this fly should be included in the model or not

%% Load data

clear all; close all;


path = uigetdir('Z:\Wilson Lab\Mel\Experiments\Uncertainty\Exp25\data\Experimental');

folderContents = dir(path);

for content = 1:length(folderContents)
   if contains(folderContents(content).name,'60D05')
       data(content) = load([folderContents(content).folder,'\',folderContents(content).name,'\analysis\summary_data.mat']);
       ver_offset(content) = load([folderContents(content).folder,'\',folderContents(content).name,'\analysis\verifying_offset.mat']);
   end
end


%% Clean and combine data

%remove empty rows
data = data(all(~cellfun(@isempty,struct2cell(data))));
ver_offset = ver_offset(~cellfun(@isempty,struct2cell(ver_offset)));

%combine data
all_offsets = [];
for fly = 1:length(data)
   all_offsets(fly,1) = data(fly).mean_reference_offset2;
   all_offsets(fly,2) = ver_offset(fly).mean_verifying_offset;
end

all_offsets_rad = deg2rad(all_offsets);

%plot both offsets
figure,
subplot(1,2,1)
plot(all_offsets,'o')
legend('ref offset','ver offset');
ylim([-180 180]);

offset_diff = rad2deg(circ_dist(all_offsets_rad(:,1),all_offsets_rad(:,2)));

subplot(1,2,2)
plot(offset_diff,'o')
ylabel('Circular difference between reference and verifier offset');
xlabel('Fly number');
ylim([-180 180]);


%% Apply criterion

%I will keep for the model those flies in which the initial and final
%offset is less than 90 deg appart.
flies_for_model = find(abs(offset_diff)<90);

%save flies for model
save('Z:\Wilson Lab\Mel\Experiments\Uncertainty\Exp25\data\Experimental\two_ND_filters_3_contrasts\flies_for_model.mat','flies_for_model')
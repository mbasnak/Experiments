%with a fly on


%load data

%DFF
meanDFF = mean(data.dff_matrix);
medianDFF = median(data.dff_matrix);

figure,
subplot(2,2,1)
plot(meanDFF)
hold on
plot(data.fr_y_ds/60 + 0.2)
title('Mean DFF across glomeruli');
ylabel('DF/F');

subplot(2,2,2)
plot(medianDFF)
hold on
plot(data.fr_y_ds/60 + 0.2)
title('Median DFF across glomeruli');


%load data

%F
for i = 1:16
    trace(i,:) = roi_data(i).trace;
end

meanTrace = mean(trace);
medianTrace = median(trace);

subplot(2,2,3)
plot(meanTrace)
hold on
plot(800000 + data.fr_y_ds*10000)
title('Mean F across glomeruli');
xlabel('Time'); ylabel('F');

subplot(2,2,4)
plot(medianTrace)
hold on
plot(600000 + data.fr_y_ds*10000)
title('Median F across glomeruli');
xlabel('Time');

%% Without a fly on

%load data (roi_data and ball data)

%F
trace = roi_data.trace;
time = linspace(0,length(trace),length(trial_time));

figure,
plot(trace)
hold on
plot(time, 3500000 + trial_bdata(:,6)*25000)
title('Mean F (smoothed)');
xlabel('Time');
ylabel('F');

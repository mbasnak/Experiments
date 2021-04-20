heading = rad2deg(data.flyPosRad');
phase = -rad2deg(data.phase);

smoothHeading = smooth(heading);
smoothPhase = smooth(phase);

figure,
subplot(5,1,1)
plot(data.time,heading,'LineWidth',1.5)
hold on
plot(data.time,phase,'LineWidth',1.5)
xlim([0, 200]);
ylim([-180 180]);
ylabel('Deg');
legend('Heading','EPG phase');

subplot(5,1,2)
plot(data.time,smoothHeading,'LineWidth',1.5)
hold on
plot(data.time,smoothPhase,'LineWidth',1.5)
xlim([0,200]);
ylim([-180 180]);
ylabel('Deg');
legend('Smooth heading','smooth EPG phase');

%Compute the offset
for frame = 1:length(heading)
    minValue(frame) = min(heading(frame),phase(frame));
    maxValue(frame) = max(heading(frame),phase(frame));
    offset(frame) = maxValue(frame) - minValue(frame);
end


subplot(5,1,3)
plot(data.time,heading,'-o','LineWidth',1.5)
hold on
plot(data.time,phase,'-o','LineWidth',1.5)
for frame = 1:length(heading)
   line([data.time(frame) data.time(frame)],[minValue(frame) maxValue(frame)],'color','k'); 
end
xlim([0,200]);
ylim([-180 180]);
ylabel('Deg');
legend('Heading','EPG phase','Offset');

subplot(5,1,4)
plot(data.time,wrapTo180(heading-phase),'-o','LineWidth',1.5,'color','k')
xlim([0,200]);
ylim([-180 180]);
ylabel('Deg');
legend('Offset');

subplot(5,1,5)
plot(data.time,smooth(wrapTo180(heading-phase)),'-o','LineWidth',1.5,'color','k')
xlim([0,200]);
ylim([-180 180]);
ylabel('Deg');
legend('Offset');

figure,
plot(data.time,smooth(wrapTo180(heading-phase)),'-o','LineWidth',1.5,'color','k')
ylim([-180 180]);
ylabel('Deg');
legend('Offset');



function plot_mvt_parameters(vel_for_deg,vel_side_deg,vel_yaw,total_mvt)

figure('Position',[100 100 1200 800]),
%forward velocity
subplot(4,5,[1 4])
plot(vel_for_deg,'color',[0.2 0.5 0.9])
title('Forward velocity');
xlim([1 length(vel_for_deg)]);

subplot(4,5,5)
histogram(vel_for_deg,'FaceColor',[0.2 0.5 0.9])

%side velocity
subplot(4,5,[6 9])
plot(vel_side_deg,'color',[0.6 0.5 0.9])
title('Side velocity');
xlim([1 length(vel_for_deg)]);

subplot(4,5,10)
histogram(vel_side_deg,'FaceColor',[0.6 0.5 0.9])

%angular velocity
subplot(4,5,[11 14])
plot(vel_yaw,'color',[0.4 0.6 0.2])
title('Angular velocity');
xlim([1 length(vel_for_deg)]);

subplot(4,5,15)
histogram(vel_yaw,'FaceColor',[0.4 0.6 0.2])

%total movement
subplot(4,5,[16 19])
plot(total_mvt,'color',[0.8 0.5 0.2])
title('Total movement');
xlim([1 length(vel_for_deg)]);

subplot(4,5,20)
histogram(total_mvt,'FaceColor',[0.8 0.5 0.2])
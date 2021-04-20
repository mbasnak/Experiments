figure,
subplot(5,1,1)
unwrapped_phase{session} = rad2deg(unwrap(-data{1,session}.data.phase));
unwrapped_heading{session} = rad2deg(unwrap(data{1,session}.data.flyPosRad));
unwrapped_stim{session} = rad2deg(unwrap(-deg2rad(data{1,session}.data.fr_ds*360/max(data{1,session}.data.fr_ds))));
plot(unwrapped_stim{session});
hold on
plot(unwrapped_heading{session});
plot(unwrapped_phase{session})
title('EPG activity phase');
legend('stimulus position','fly heading','phase position');

subplot(5,1,2)
%calculate velocity roughly per second
relevant_timepoints = 10:9:276;
for timepoint = 10:9:276
    stim_vel(timepoint) = unwrapped_stim{session}(timepoint)-unwrapped_stim{session}(timepoint-9);
    fly_vel(timepoint) = unwrapped_heading{session}(timepoint)-unwrapped_heading{session}(timepoint-9);
    phase_vel(timepoint) = unwrapped_phase{session}(timepoint)-unwrapped_phase{session}(timepoint-9);
end

stim_vel = stim_vel(relevant_timepoints);
fly_vel = fly_vel(relevant_timepoints);
phase_vel = phase_vel(relevant_timepoints);
plot(stim_vel)
hold on
plot(fly_vel)
plot(phase_vel)

%linearly interpolate to obtain more points
subplot(5,1,3)
stim_vel_interp = interp1(1:length(stim_vel),stim_vel,1:length(unwrapped_stim{session}));
fly_vel_interp = interp1(1:length(fly_vel),fly_vel,1:length(unwrapped_heading{session}));
phase_vel_interp = interp1(1:length(phase_vel),phase_vel,1:length(unwrapped_phase{session}));
plot(stim_vel_interp)
hold on
plot(fly_vel_interp)
plot(phase_vel_interp)

%smooth data
subplot(5,1,4)
smoothed_stim_vel = smooth(stim_vel_interp,'rlowess');
smoothed_fly_vel = fly_vel_interp;%smooth(fly_vel_interp,'rloess');
smoothed_phase_vel = phase_vel_interp;%;smoothdata(phase_vel_interp);

plot(smoothed_stim_vel)
hold on
plot(smoothed_fly_vel)
plot(smoothed_phase_vel)

subplot(5,1,5)
plot(abs(smoothed_stim_vel))
hold on
plot(abs(smoothed_fly_vel))
plot(abs(smoothed_phase_vel))



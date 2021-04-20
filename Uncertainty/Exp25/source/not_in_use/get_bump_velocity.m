%% Take good tracking portion to extract velocity smoothing parameters

%import data

load('Z:\Wilson Lab\Mel\Experiments\Uncertainty\Exp25\data\Experimental\two_ND_filters_3_contrasts\20201019_60D05_7f\analysis\analysis_sid_1_tid_0.mat');
%Obtain the frames where the stim changes using the derivative of the
%signal from the y dimension of the panels
changeContrast = find(abs(diff(data.fr_y_ds))>1);

%% Identify contrast order

%Check the position function
pos_function = data.run_obj.function_number;

%Define the order of intensities according to the function used
if pos_function == 195  
    Intensities = [1,2,1,3,2,3];
elseif pos_function == 196  
    Intensities = [2,1,3,1,2,3];
else
    Intensities = [3,1,2,1,2,3];
end
%% Set block limits

blockLimits{1} = [1,changeContrast(1)-1];

if length(changeContrast) == 5
    for block = 2:5
        blockLimits{block} = [changeContrast(block-1),changeContrast(block)-1];
    end
    blockLimits{6} = [changeContrast(5),length(heading)];
elseif length(changeContrast) == 6
    for block = 2:6
        blockLimits{block} = [changeContrast(block-1),changeContrast(block)-1];
    end
else
    for block = 2:4
        blockLimits{block} = [changeContrast(block-1),changeContrast(block)-1];
    end
     blockLimits{5} = [changeContrast(4),length(heading)];
end

%Find limits of last high contrast bout
high_contrast_bouts = find(Intensities==3);
last_high_contrast = high_contrast_bouts(end);

%% Get and plot bouts of good tracking

example_phase = wrapTo360(rad2deg(data.dff_pva(blockLimits{high_contrast_bouts(1)}(1)+500:blockLimits{high_contrast_bouts(1)}(2))));
example_stim = data.panel_angle(blockLimits{high_contrast_bouts(1)}(1)+500:blockLimits{high_contrast_bouts(1)}(2));
example_offset = data.offset(blockLimits{high_contrast_bouts(1)}(1)+500:blockLimits{high_contrast_bouts(1)}(2));

figure,
subplot(2,1,1)
plot(example_phase)
hold on
plot(example_stim)
subplot(2,1,2)
plot(example_offset)

%% Obtain velocity

%get angular velocity
figure('Position',[100 100 1600 1000]),
%plot position
subplot(5,4,[1 3])
plot(example_phase)
hold on
plot(example_stim)
legend('EPG phase','Stimulus position')
title('Bump and stimulus position');
xlim([0 length(example_phase)]);

%unwrap position
subplot(5,4,[5 7])
unwrapped_example_phase = unwrap(deg2rad(example_phase));
unwrapped_example_stim = unwrap(deg2rad(example_stim));
plot(unwrapped_example_phase)
hold on
plot(unwrapped_example_stim)
legend('unwrapped EPG phase','unwrapped stimulus position')
title('Unwrapped bump and stimulus position');
xlim([0 length(example_phase)]);

%smooth position
subplot(5,4,[9 11])
smooth_example_phase = smooth(unwrapped_example_phase);
smooth_example_stim = smooth(unwrapped_example_stim);
plot(smooth_example_phase)
hold on
plot(smooth_example_stim)
legend('smooth unwrapped EPG phase','smooth unwrapped stimulus position')
title('Smooth unwrapped bump and stimulus position');
xlim([0 length(example_phase)]);

%Plot positional correlation
subplot(5,4,[4,8,12])
scatter(smooth_example_phase,smooth_example_stim)
title('Correlation between unwrapped bump and stimulus position');

%take derivative
subplot(5,4,[13 15])
diff_example_phase = diff(smooth_example_phase);
diff_example_stim = diff(smooth_example_stim);
plot(diff_example_phase)
hold on
plot(diff_example_stim)
legend('EPG phase vel','Stimulus vel')
title('Bump and stimulus velocity');
xlim([0 length(example_phase)]);

%smooth velocity
subplot(5,4,[17 19])
smooth_diff_example_phase = smooth(diff_example_phase);
smooth_diff_example_stim = smooth(diff_example_stim);
plot(smooth_diff_example_phase)
hold on
plot(smooth_diff_example_stim)
legend('smooth EPG phase vel','smooth stimulus vel')
title('Smooth bump and stimulus velocity');
xlim([0 length(example_phase)]);

subplot(5,4,[16,20])
scatter(smooth_diff_example_phase,smooth_diff_example_stim)
title('Correlation between bump and stimulus velocity');

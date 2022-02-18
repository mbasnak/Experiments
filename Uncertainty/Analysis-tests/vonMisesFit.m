%Code to compare different options for fitting the PB data

clear all; close all;
%% Load sample data

load('continuous_data.mat')

%% 1) Function we had been using so far

%1) Convert the distances along the PB midline to angular distances, taking
%into account that the full length of the PB corresponds to 2pi + (7/8)pi
%in angular distances (since we are grouping the middle glomeruli)
angular_midline_distances = midline_distances*(2*pi + 2*pi*(7/8))/max(midline_distances);
%put back on circle
angular_midline_distances_2pi = mod(angular_midline_distances,2*pi);

%2) Create the fit function
ft = fittype('a*exp(k*cos(x-u))');

for timepoint = 1:length(dff)
    %3) Fit the data
    % add the minimum value
    offset = min(dff(timepoint,:));
    data_to_fit = dff(timepoint,:) - offset;
    [model_data,gof] = fit(angular_midline_distances_2pi,data_to_fit',ft,'StartPoint',[2,2,4]);
    adj_rs(timepoint) = gof.adjrsquare;
    
    %4) Get the coefficient values.
    coeff_names = coeffnames(ft);
    coefficientValues = coeffvalues(model_data);
    a = coefficientValues(strcmp(coeff_names,'a'));
    k = coefficientValues(strcmp(coeff_names,'k'));
    u = coefficientValues(strcmp(coeff_names,'u'));
    %5) Correct for negative k values (ambiguity in fit)
    if k < 0
        k = -k;
        u = u + pi;
    end
    bump_pos(timepoint) = mod(u,2*pi);
    
    %6) Get bump parameters
    bump_mag(timepoint) = a * exp(k);
    bump_width(timepoint) = 2 * abs(acos(1-log(2)/k));
    
    
    %7) Plot the original data and the fit
%     figure,
%     subplot(2,1,1)
%     plot(angular_midline_distances,dff_data(timepoint,:)')
%     hold on
%     plot(angular_midline_distances,feval(model_data,angular_midline_distances)+offset)
%     %add the bump position estimate
%     plot(u,feval(model_data,u)+offset,'ro')
%     xlabel('Angular distance (radians)');
%     ylabel('DF/F');
%     title(['Frame #',num2str(timepoint),' AdjR2 =',num2str(gof.adjrsquare)]);
%     legend('data','fit');
%     %Add the bump parameters
%     subplot(2,1,2)
%     text(0,0.5,['Bump magnitude = ',num2str(bump_mag(timepoint))]);
%     hold on
%     text(0.5,0.5,['Bump width = ',num2str(bump_width(timepoint))]);
%     text(0.25,0,['Bump pos = ',num2str(u)]); axis off
    
end

%% 2) If I don't correct for negative k values

for timepoint = 1:length(dff)
    %3) Fit the data
    % add the minimum value
    offset = min(dff(timepoint,:));
    data_to_fit = dff(timepoint,:) - offset;
    [model_data,gof] = fit(angular_midline_distances_2pi,data_to_fit',ft,'StartPoint',[2,2,4]);
    adj_rs(timepoint) = gof.adjrsquare;
    
    %4) Get the coefficient values.
    coeff_names = coeffnames(ft);
    coefficientValues = coeffvalues(model_data);
    a = coefficientValues(strcmp(coeff_names,'a'));
    k = coefficientValues(strcmp(coeff_names,'k'));
    u = coefficientValues(strcmp(coeff_names,'u'));

    uncorrected_bump_pos(timepoint) = mod(u,2*pi);
    
    %5) Get bump parameters
    uncorrected_bump_mag(timepoint) = a * exp(k);
    uncorrected_bump_width(timepoint) = 2 * abs(acos(1-log(2)/k)); 
   
end

%% Compare the bump parameters in both cases

figure,
subplot(3,1,1)
plot(bump_mag)
hold on
plot(uncorrected_bump_mag)
legend('Corrected for -k','Uncorrected');
title('Bump magnitude');

subplot(3,1,2)
plot(bump_width)
hold on
plot(uncorrected_bump_width)
legend('Corrected for -k','Uncorrected');
title('Bump width');

subplot(3,1,3)
plot(bump_pos)
hold on
plot(uncorrected_bump_pos)
legend('Corrected for -k','Uncorrected');
title('Bump position');

saveas(gcf,'Z:\Wilson Lab\Mel\Experiments\Uncertainty\Analysis-tests\correcting_for_k.png')
%The three plots differ in several places


%% 3) Using the 'correct' equation, with takes into account the offset using a constant c

%1) Create the fit function
ft_with_constant = fittype('a*exp(k*cos(x-u))+c');

for timepoint = 1:length(dff)
    %2) Fit the data
    data_to_fit = dff(timepoint,:);
    [model_data,gof] = fit(angular_midline_distances_2pi,data_to_fit',ft_with_constant,'StartPoint',[2,2,4,0]);
    adj_rs(timepoint) = gof.adjrsquare;
    
    %3) Get the coefficient values.
    coeff_names = coeffnames(ft);
    coefficientValues = coeffvalues(model_data);
    a = coefficientValues(strcmp(coeff_names,'a'));
    k = coefficientValues(strcmp(coeff_names,'k'));
    u = coefficientValues(strcmp(coeff_names,'u'));
    c = coefficientValues(strcmp(coeff_names,'c'));
    
    %4) Correct for negative k values (ambiguity in fit)
    if k < 0
        k = -k;
        u = u + pi;
    end
    bump_pos_with_constant(timepoint) = mod(u,2*pi);
    
    %5) Get bump parameters
    bump_mag_with_constant(timepoint) = a * (exp(k) - exp(-k));
    %bump_width_with_constant(timepoint) = 2 * abs(acos(1-log(2)/k));
       
end

%% Compare the constant method with the original one we used

figure,
subplot(3,1,1)
plot(bump_mag)
hold on
plot(bump_mag_with_constant)
legend('Without constant','With constant');
title('Bump magnitude');

subplot(3,1,2)
plot(bump_width)
% hold on
% plot(bump_width_with_constant)
%legend('Without constant','With constant');
title('Bump width');

subplot(3,1,3)
plot(bump_pos)
hold on
plot(bump_pos_with_constant)
legend('Without constant','With constant');
title('Bump position');

saveas(gcf,'Z:\Wilson Lab\Mel\Experiments\Uncertainty\Analysis-tests\with_and_without_constant.png')

%% 4) Using a von Mises profile function (centered around zero) to obtain the bump parameters

%1) Create the fit function
ft_with_constant = fittype('a*exp(k*cos(x-u))+c');

for timepoint = 1:length(dff)
    %2) Fit the data
    data_to_fit = dff(timepoint,:);
    [model_data,gof] = fit(angular_midline_distances_2pi,data_to_fit',ft_with_constant,'StartPoint',[2,2,4,0]);
    adj_rs(timepoint) = gof.adjrsquare;
    
    %3) Get the coefficient values.
    coeff_names = coeffnames(ft);
    coefficientValues = coeffvalues(model_data);
    a = coefficientValues(strcmp(coeff_names,'a'));
    k = coefficientValues(strcmp(coeff_names,'k'));
    u = coefficientValues(strcmp(coeff_names,'u'));
    c = coefficientValues(strcmp(coeff_names,'c'));
    
    %4) Correct for negative k values (ambiguity in fit)
    if k < 0
        k = -k;
        u = u + pi;
    end
    bump_pos_with_constant(timepoint) = mod(u,2*pi);
    
    %5) Plug the coefficients into a von Mises pdf
    %[p, alpha] = circ_vmpdf(angular_midline_distances_2pi, u, k);
    [p, alpha] = circ_vmpdf([0:0.006:2*pi], u, k);    
    
    % 5) rescale von mises function to match orignial data
    pMinSubtract = p - min(p); % min subtract
    p_norm = pMinSubtract / max( pMinSubtract ); % max divide
    p_rescaled = ( p_norm *  max(data_to_fit) ) + min(data_to_fit);

    %6) Compute bump parameters from that fit
    %Bump magnitude
    von_mises_bump_mag(timepoint) = max(p_rescaled)-min(p_rescaled);
    
    %Bump width
    %Find the half max value in the y axis as the middle point between the min
    %and max y axis values
    half_max = (max(p_rescaled)-min(p_rescaled))/2 + min(p_rescaled);
    [ex_bump_mag_interp I_interp] = max(p_rescaled);
    %Find in each half the index closest to the half max value
    diff_data = abs(p_rescaled-half_max);
    [sortedVals,indexes] = sort(diff_data);
    diff_indexes = abs(indexes-indexes(1));
    %Remove all the indexes that are less than 175 datapoints away (i.e., a little more than 1 glomerulus) from
    %index(1)
    indexes(diff_indexes<195 & diff_indexes>0)=NaN;
    indexes = indexes(~isnan(indexes));
    two_indexes = [indexes(1), indexes(2)];
    I1 = min(two_indexes);
    I2 = max(two_indexes);
    if (all(two_indexes>I_interp) | all(two_indexes<I_interp))
        half_max_w = I1+length(p_rescaled)-I2;
    else
        half_max_w = I2-I1;
    end
    %Convert to EB coordinates
    von_mises_bump_width(timepoint) = half_max_w*8/length(p_rescaled);
     
end

%% Compare bump parameters with original method


figure,
subplot(2,1,1)
plot(bump_mag)
hold on
plot(von_mises_bump_mag)
legend('Original method','with von Mises');
title('Bump magnitude');

subplot(2,1,2)
plot(bump_width)
hold on
plot(von_mises_bump_width)
legend('Original method','with von Mises');
title('Bump width');

saveas(gcf,'Z:\Wilson Lab\Mel\Experiments\Uncertainty\Analysis-tests\from_von_Mises.png')

%% Example timepoint

figure,
subplot(2,1,1)
plot(data_to_fit)
title('DF/F data')

subplot(2,1,2)
plot(p_rescaled)
hold on
plot(I1,p_rescaled(I1),'ro')
plot(I2,p_rescaled(I2),'ro')
line([I1 I2],[p_rescaled(I1) p_rescaled(I1)],'color','r')
title('Bump width')

saveas(gcf,'Z:\Wilson Lab\Mel\Experiments\Uncertainty\Analysis-tests\example_timepoint.png')



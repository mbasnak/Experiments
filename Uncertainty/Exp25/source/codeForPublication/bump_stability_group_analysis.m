%Code to analyze the bump stability in the closed-loop bouts, and compare
%it to that of a laser power control

clear all; close all;

%% Load data

path = uigetdir('Z:\Wilson Lab\Mel\Experiments\Uncertainty\Exp25\data\Experimental');

folderContents = dir(path);

for content = 1:length(folderContents)
   if contains(folderContents(content).name,'60D05')
       data(content) = load([folderContents(content).folder,'\',folderContents(content).name,'\analysis\summary_data.mat']);
   end
end

%Remove empty rows
data = data(all(~cellfun(@isempty,struct2cell(data))));

%% Bump stability analysis using a movement threshold

%Define the different bin numbers
bins = [5,10,15,20,30,40,50,100];
% 
% all_bm_darkness = [];
% all_hw_darkness = [];
% all_bm_lc = [];
% all_hw_lc = [];
% all_bm_hc = [];
% all_hw_hc = [];
% 
% figure('Position',[100, 100, 1600, 1000]),
% for fly = 1:length(data)
%     subplot(3,2,1)
%     for bin = 1:length(bins)
%         if data(fly).pval_bm_thresh{1,1}(bin) < 0.05
%             plot(bins(bin),data(fly).correlations_bm_thresh{1,1}(bin),'ro','MarkerSize',4)
%         else
%             plot(bins(bin),data(fly).correlations_bm_thresh{1,1}(bin),'ko','MarkerSize',4)
%         end
%         hold on
%     end
%     plot(bins,data(fly).correlations_bm_thresh{1,1},'color',[.5 .5 .5])
%     ylim([-1 1]);
%     title('Relationship between bump magnitude and bump instability (bump var)');
%     ylabel('Darkness bouts','Rotation',0,'fontweight','bold');
%     xh = get(gca,'ylabel'); % handle to the label object
%     p = get(xh,'position'); % get the current position property
%     p(1) = 1.1*p(1) ;        % increase the distance,
%     set(xh,'position',p);   % set the new position
%     %add mean
%     all_bm_darkness = [all_bm_darkness;data(fly).correlations_bm_thresh{1,1}];
%     if fly == length(data)
%         plot(bins,mean(all_bm_darkness),'-ko','markerfacecolor','k','linewidth',3)
%     end
%     %add zero line
%     line([0 100],[0 0],'color','b','linewidth',2);
%     
%     subplot(3,2,2)
%     for bin = 1:length(bins)
%         if data(fly).pval_hw_thresh{1,1}(bin) < 0.05
%             plot(bins(bin),data(fly).correlations_hw_thresh{1,1}(bin),'ro','MarkerSize',4)
%         else
%             plot(bins(bin),data(fly).correlations_hw_thresh{1,1}(bin),'ko','MarkerSize',4)
%         end
%         hold on
%     end
%     plot(bins,data(fly).correlations_hw_thresh{1,1},'color',[.5 .5 .5])
%     ylim([-1 1]);
%     title('Relationship between bump width and bump instability (bump var)');
%     %add mean
%     all_hw_darkness = [all_hw_darkness;data(fly).correlations_hw_thresh{1,1}];
%     if fly == length(data)
%         plot(bins,mean(all_hw_darkness),'-ko','markerfacecolor','k','linewidth',3)
%     end
%     %add zero line
%     line([0 100],[0 0],'color','b','linewidth',2);
%     
%     subplot(3,2,3)
%     for bin = 1:length(bins)
%         if data(fly).pval_bm_thresh{1,2}(bin) < 0.05
%             plot(bins(bin),data(fly).correlations_bm_thresh{1,2}(bin),'ro','MarkerSize',4)
%         else
%             plot(bins(bin),data(fly).correlations_bm_thresh{1,2}(bin),'ko','MarkerSize',4)
%         end
%         hold on
%     end
%     plot(bins,data(fly).correlations_bm_thresh{1,2},'color',[.5 .5 .5])
%     ylim([-1 1]);
%     ylabel('Low contrast bouts','Rotation',0,'fontweight','bold');
%     xh = get(gca,'ylabel'); % handle to the label object
%     p = get(xh,'position'); % get the current position property
%     p(1) = 1.1*p(1) ;        % increase the distance,
%     set(xh,'position',p);   % set the new position
%     %add mean
%     all_bm_lc = [all_bm_lc;data(fly).correlations_bm_thresh{1,2}];
%     if fly == length(data)
%         plot(bins,mean(all_bm_lc),'-ko','markerfacecolor','k','linewidth',3)
%     end
%     %add zero line
%     line([0 100],[0 0],'color','b','linewidth',2);
%     
%     subplot(3,2,4)
%     for bin = 1:length(bins)
%         if data(fly).pval_hw_thresh{1,2}(bin) < 0.05
%             plot(bins(bin),data(fly).correlations_hw_thresh{1,2}(bin),'ro','MarkerSize',4)
%         else
%             plot(bins(bin),data(fly).correlations_hw_thresh{1,2}(bin),'ko','MarkerSize',4)
%         end
%         hold on
%     end
%     plot(bins,data(fly).correlations_hw_thresh{1,2},'color',[.5 .5 .5])
%     ylim([-1 1]);
%     %add mean
%     all_hw_lc = [all_hw_lc;data(fly).correlations_hw_thresh{1,2}];
%     if fly == length(data)
%         plot(bins,mean(all_hw_lc),'-ko','markerfacecolor','k','linewidth',3)
%     end
%     %add zero line
%     line([0 100],[0 0],'color','b','linewidth',2);
%     
%     subplot(3,2,5)
%     for bin = 1:length(bins)
%         if data(fly).pval_bm_thresh{1,3}(bin) < 0.05
%             plot(bins(bin),data(fly).correlations_bm_thresh{1,3}(bin),'ro','MarkerSize',4)
%         else
%             plot(bins(bin),data(fly).correlations_bm_thresh{1,3}(bin),'ko','MarkerSize',4)
%         end
%         hold on
%     end
%     plot(bins,data(fly).correlations_bm_thresh{1,3},'color',[.5 .5 .5])
%     ylim([-1 1]);
%     xlabel('Bin number');
%     ylabel('High contrast bouts','Rotation',0,'fontweight','bold');
%     xh = get(gca,'ylabel'); % handle to the label object
%     p = get(xh,'position'); % get the current position property
%     p(1) = 1.1*p(1) ;        % increase the distance,
%     set(xh,'position',p);   % set the new position
%     %add mean
%     all_bm_hc = [all_bm_hc;data(fly).correlations_bm_thresh{1,3}];
%     if fly == length(data)
%         plot(bins,mean(all_bm_hc),'-ko','markerfacecolor','k','linewidth',3)
%     end
%     %add zero line
%     line([0 100],[0 0],'color','b','linewidth',2);
%     
%     subplot(3,2,6)
%     for bin = 1:length(bins)
%         if data(fly).pval_hw_thresh{1,3}(bin) < 0.05
%             plot(bins(bin),data(fly).correlations_hw_thresh{1,3}(bin),'ro','MarkerSize',4)
%         else
%             plot(bins(bin),data(fly).correlations_hw_thresh{1,3}(bin),'ko','MarkerSize',4)
%         end
%         hold on
%     end
%     plot(bins,data(fly).correlations_hw_thresh{1,3},'color',[.5 .5 .5])
%     ylim([-1 1]);
%     xlabel('Bin number');
%     %add mean
%     all_hw_hc = [all_hw_hc;data(fly).correlations_hw_thresh{1,3}];
%     if fly == length(data)
%         plot(bins,mean(all_hw_hc),'-ko','markerfacecolor','k','linewidth',3)
%     end
%     %add zero line
%     line([0 100],[0 0],'color','b','linewidth',2);
%     
% end
% 
% %save figure
% saveas(gcf,[path,'\globalPlots\bumpStabilityThresh_betweenBlocks.png']);
% 
% 
% %% Model for bump stability
% 
% for bin = 1:length(bins)
%     
%     for contrast = 1:3
%         %create tables with the relevant variables
%         bsData_per_contrast{contrast,bin} = array2table(zeros(0,3),'VariableNames', {'bump_var','bump_mag','bump_width'});
%         Fly = [];
%         for fly = 1:length(data)
%             Fly = [Fly,repelem(fly,length(data(fly).bump_var_per_contrast{contrast,bin}))];
%             bs_per_contrast = array2table([data(fly).bump_var_per_contrast{contrast,bin}',data(fly).binned_bump_mag_per_contrast{contrast,bin}',data(fly).binned_half_width_per_contrast{contrast,bin}'],'VariableNames', {'bump_var','bump_mag','bump_width'});
%             bsData_per_contrast{contrast,bin} = [bsData_per_contrast{contrast,bin};bs_per_contrast];
%         end
%         bsData_per_contrast{contrast,bin} = addvars(bsData_per_contrast{contrast,bin},Fly');
%         bsData_per_contrast{contrast,bin}.Properties.VariableNames{'Var4'} = 'Fly';
%         
%         mdl_bs_bm_per_contrast{contrast,bin} = fitlme(bsData_per_contrast{contrast,bin},'bump_var~1+bump_mag+(1|Fly)');
%         coef_bm_per_contrast(contrast,bin) = mdl_bs_bm_per_contrast{contrast,bin}.Coefficients.Estimate(2);
%         pval_bm_per_contrast(contrast,bin) = mdl_bs_bm_per_contrast{contrast,bin}.Coefficients.pValue(2);
%         
%         mdl_bs_hw_per_contrast{contrast,bin} = fitlme(bsData_per_contrast{contrast,bin},'bump_var~1+bump_width+(1|Fly)');
%         coef_hw_per_contrast(contrast,bin) = mdl_bs_hw_per_contrast{contrast,bin}.Coefficients.Estimate(2);
%         pval_hw_per_contrast(contrast,bin) = mdl_bs_hw_per_contrast{contrast,bin}.Coefficients.pValue(2);
%     end
%     
% end
% 
% %repeat figure with these coefficients
% figure('Position',[100,100,1400,1000]),
% for contrast = 1:3
%     subplot(3,2,contrast*2-1)
%     plot(bins,coef_bm_per_contrast(contrast,:),'color',[.5 .5 .5],'linewidth',2)
%     hold on
%     for bin = 1:length(bins)
%         if pval_bm_per_contrast(contrast,bin) < 0.05
%             plot(bins(bin),coef_bm_per_contrast(contrast,bin),'ro','linewidth',2)
%         else
%             plot(bins(bin),coef_bm_per_contrast(contrast,bin),'ko','linewidth',2)
%         end
%         hold on
%     end
%     ylim([-1 1]);
%     title('Relationship between bump magnitude and bump instability (bump var)');
%     ylabel('Regression coefficient');
%     xlabel('Bin number');
%     line([0 100],[0 0],'color','b','linewidth',2);
%     
%     subplot(3,2,contrast*2)
%     plot(bins,coef_hw_per_contrast(contrast,:),'color',[.5 .5 .5],'linewidth',2)
%     hold on
%     for bin = 1:length(bins)
%         if pval_hw_per_contrast(contrast,bin) < 0.05
%             plot(bins(bin),coef_hw_per_contrast(contrast,bin),'ro','linewidth',2)
%         else
%             plot(bins(bin),coef_hw_per_contrast(contrast,bin),'ko','linewidth',2)
%         end
%         hold on
%     end
%     ylim([-1 1]);
%     title('Relationship between bump width and bump instability (bump var)');
%     ylabel('Regression coefficient');
%     xlabel('Bin number');
%     line([0 100],[0 0],'color','b','linewidth',2);
%     
%   
% end
% 
% %save figure
% saveas(gcf,[path,'\globalPlots\bumpStabilityThresh_with_intercept_pval.png']);
% % 
% % %% Model for bump stability without intercept
% % 
% % for bin = 1:length(bins)
% %     
% %     for contrast = 1:3
% %         %create tables with the relevant variables
% %         bsData_per_contrast{contrast,bin} = array2table(zeros(0,3),'VariableNames', {'bump_var','bump_mag','bump_width'});
% %         Fly = [];
% %         for fly = 1:length(data)
% %             Fly = [Fly,repelem(fly,length(data(fly).bump_var_per_contrast{contrast,bin}))];
% %             bs_per_contrast = array2table([data(fly).bump_var_per_contrast{contrast,bin}',data(fly).binned_bump_mag_per_contrast{contrast,bin}',data(fly).binned_half_width_per_contrast{contrast,bin}'],'VariableNames', {'bump_var','bump_mag','bump_width'});
% %             bsData_per_contrast{contrast,bin} = [bsData_per_contrast{contrast,bin};bs_per_contrast];
% %         end
% %         bsData_per_contrast{contrast,bin} = addvars(bsData_per_contrast{contrast,bin},Fly');
% %         bsData_per_contrast{contrast,bin}.Properties.VariableNames{'Var4'} = 'Fly';
% %         
% %         mdl_bs_bm2{contrast,bin} = fitlme(bsData_per_contrast{contrast,bin},'bump_var~-1+bump_mag+(1|Fly)');
% %         coef_bm_per_contrast(contrast,bin) = mdl_bs_bm2{contrast,bin}.Coefficients.Estimate;
% %         pval_bm_per_contrast(contrast,bin) = mdl_bs_bm2{contrast,bin}.Coefficients.pValue;
% %         
% %         mdl_bs_hw2{contrast,bin} = fitlme(bsData_per_contrast{contrast,bin},'bump_var~-1+bump_width+(1|Fly)');
% %         coef_hw_per_contrast(contrast,bin) = mdl_bs_hw2{contrast,bin}.Coefficients.Estimate;
% %         pval_hw_per_contrast(contrast,bin) = mdl_bs_hw2{contrast,bin}.Coefficients.pValue;
% %     end
% %     
% % end
% % 
% % %repeat figure with these coefficients
% % figure('Position',[100,100,1400,1000]),
% % for contrast = 1:3
% %     subplot(3,2,contrast*2-1)
% %     plot(bins,coef_bm_per_contrast(contrast,:),'color',[.5 .5 .5],'linewidth',2)
% %     hold on
% %     for bin = 1:length(bins)
% %         if pval_bm_per_contrast(contrast,bin) < 0.05
% %             plot(bins(bin),coef_bm_per_contrast(contrast,bin),'ro','linewidth',2)
% %         else
% %             plot(bins(bin),coef_bm_per_contrast(contrast,bin),'ko','linewidth',2)
% %         end
% %         hold on
% %     end
% %     ylim([-1 1]);
% %     title('Relationship between bump magnitude and bump instability (bump var)');
% %     ylabel('Regression coefficient');
% %     xlabel('Bin number');
% %     line([0 100],[0 0],'color','b','linewidth',2);
% %     
% %     subplot(3,2,contrast*2)
% %     plot(bins,coef_hw_per_contrast(contrast,:),'color',[.5 .5 .5],'linewidth',2)
% %     hold on
% %     for bin = 1:length(bins)
% %         if pval_hw_per_contrast(contrast,bin) < 0.05
% %             plot(bins(bin),coef_hw_per_contrast(contrast,bin),'ro','linewidth',2)
% %         else
% %             plot(bins(bin),coef_hw_per_contrast(contrast,bin),'ko','linewidth',2)
% %         end
% %         hold on
% %     end
% %     ylim([-1 1]);
% %     title('Relationship between bump width and bump instability (bump var)');
% %     ylabel('Regression coefficient');
% %     xlabel('Bin number');
% %     line([0 100],[0 0],'color','b','linewidth',2);
% %     
% %   
% % end
% % 
% % %save figure
% % saveas(gcf,[path,'\globalPlots\bumpStabilityThresh_without_intercept_pval.png']);
% % 
% % 
% % %% Compare both models
% % 
% % %we're going to compare the models with and without intercept looking at
% % %deviance and AIC
% % for bin = 1:length(bins)
% %     for contrast = 1:3
% %     deviance_bm(bin,contrast,1) = mdl_bs_bm{contrast,bin}.ModelCriterion.Deviance; %bm, with intercept
% %     deviance_bm(bin,contrast,2) = mdl_bs_bm2{contrast,bin}.ModelCriterion.Deviance; %bm, without intercept
% %     deviance_hw(bin,contrast,1) = mdl_bs_hw{contrast,bin}.ModelCriterion.Deviance; %hw, with intercept
% %     deviance_hw(bin,contrast,2) = mdl_bs_hw2{contrast,bin}.ModelCriterion.Deviance; %hw, without intercept
% %     end
% % end
% % 
% % figure('Position',[100 100 1400 1000]),
% % for contrast = 1:3
% % subplot(3,2,contrast*2-1)
% % plot(bins,deviance_bm(:,contrast,1),'-o')
% % hold on
% % plot(bins,deviance_bm(:,contrast,2),'-o')
% % title('Deviance comparison for bump mag models');
% % legend('With intercept','Without intercept');
% % 
% % 
% % subplot(3,2,contrast*2)
% % plot(bins,deviance_hw(:,contrast,1),'-o')
% % hold on
% % plot(bins,deviance_hw(:,contrast,2),'-o')
% % title('Deviance comparison for bump width models');
% % legend('With intercept','Without intercept');
% % 
% % end
% % 
% % %save figure
% % saveas(gcf,[path,'\globalPlots\modelComparison_between_blocks.png']);
% 
% 
% %% Look at the magnitude of the effect on bump variation for the between block comparison
% 
% %look at the estimates in the model and how much they're affected by
% %changes in bump magnitude or bump width
% for bin = 1:length(bins)
%     for contrast = 1:3
%         bump_var_BM_one_per_contrast(contrast,bin) = mdl_bs_bm_per_contrast{contrast,bin}.Coefficients.Estimate(1)+mdl_bs_bm_per_contrast{contrast,bin}.Coefficients.Estimate(2);
%         bump_var_BM_two_per_contrast(contrast,bin) = mdl_bs_bm_per_contrast{contrast,bin}.Coefficients.Estimate(1)+2*mdl_bs_bm_per_contrast{contrast,bin}.Coefficients.Estimate(2);
%         
%         bump_var_HW_two_per_contrast(contrast,bin) = mdl_bs_hw_per_contrast{contrast,bin}.Coefficients.Estimate(1)+2*mdl_bs_hw_per_contrast{contrast,bin}.Coefficients.Estimate(2);
%         bump_var_HW_three_per_contrast(contrast,bin) = mdl_bs_hw_per_contrast{contrast,bin}.Coefficients.Estimate(1)+3*mdl_bs_hw_per_contrast{contrast,bin}.Coefficients.Estimate(2);
%     end
% end
% 
% figure('Position',[100 100 1600 1000]),
% for contrast = 1:3
%     %Plot the percentage change
%     subplot(3,2,contrast*2-1)
%     change_BM = 100*(bump_var_BM_two_per_contrast(contrast,:)-bump_var_BM_one_per_contrast(contrast,:))./bump_var_BM_one_per_contrast(contrast,:);
%     plot(bins,change_BM,'-ko');
%     hold on
%     line([0 100],[0 0],'linewidth',2,'color','b');
%     xlabel('Bin number');
%     ylabel('Percentage change in bump var');
%     title('Change in bump var per unit change in BM');
%     
%     subplot(3,2,contrast*2)
%     change_HW = 100*(bump_var_HW_three_per_contrast(contrast,:)-bump_var_HW_two_per_contrast(contrast,:))./bump_var_HW_two_per_contrast(contrast,:);
%     plot(bins,change_HW,'-ko');
%     hold on
%     line([0 100],[0 0],'linewidth',2,'color','b');
%     title('Change in bump var per unit change in HW');
%     xlabel('Bin number');
% end
% 
% 
% %save figure
% saveas(gcf,[path,'\globalPlots\bumpStability_betweenBlocks_effect_magnitude.png']);
% 
% %% Bump stability vs bump parameters across stimulus types
% 
% all_bm = [];
% all_hw = [];
% 
% figure('Position',[100,100,1400,1000]),
% for fly = 1:length(data)
%     subplot(1,2,1)
%     for bin = 1:length(bins)
%         if data(fly).pval_bm_thresh_pooled(bin) < 0.05
%             plot(bins(bin),data(fly).correlations_bm_thresh_pooled(bin),'ro','MarkerSize',4)
%         else
%             plot(bins(bin),data(fly).correlations_bm_thresh_pooled(bin),'ko','MarkerSize',4)
%         end
%         hold on
%     end
%     plot(bins,data(fly).correlations_bm_thresh_pooled,'color',[.5 .5 .5])
%     %add mean
%     all_bm = [all_bm;data(fly).correlations_bm_thresh_pooled];
%     if fly == length(data)
%         plot(bins,mean(all_bm),'-ko','markerfacecolor','k','linewidth',3)
%     end
%     ylim([-1 1]);
%     title('Relationship between bump magnitude and bump instability (bump var)');
%     ylabel('Correlation coefficient');
%     xlabel('Bin number');
%     line([0 100],[0 0],'color','b','linewidth',2);
%     
%     subplot(1,2,2)
%     for bin = 1:length(bins)
%         if data(fly).pval_hw_thresh_pooled(bin) < 0.05
%             plot(bins(bin),data(fly).correlations_hw_thresh_pooled(bin),'ro','MarkerSize',4)
%         else
%             plot(bins(bin),data(fly).correlations_hw_thresh_pooled(bin),'ko','MarkerSize',4)
%         end
%         hold on
%     end
%     plot(bins,data(fly).correlations_hw_thresh_pooled,'color',[.5 .5 .5])
%     %add mean
%     all_hw = [all_hw;data(fly).correlations_hw_thresh_pooled];
%     if fly == length(data)
%         plot(bins,mean(all_hw),'-ko','markerfacecolor','k','linewidth',3)
%     end
%     ylim([-1 1]);
%     title('Relationship between bump width and bump instability (bump var)');
%     ylabel('Correlation coefficient');
%     xlabel('Bin number');
%     line([0 100],[0 0],'color','b','linewidth',2);
% end
% 
% %save figure
% saveas(gcf,[path,'\globalPlots\bumpStabilityThresh_acrossBlocks.png']);
% 
% %% Run models of bump stability for different bin numbers
% 
% for bin = 1:length(bins)
%     
%     %combine the tables into one
%     bsData{bin} = array2table(zeros(0,3),'VariableNames', {'bump_var','bump_mag','bump_width'});
%     Fly = [];
%     for fly = 1:length(data)
%         Fly = [Fly,repelem(fly,length(data(fly).bs_table{1,bin}.bump_mag))];
%         bsData{bin} = [bsData{bin};data(fly).bs_table{1,bin}];
%     end
%     bsData{bin} = addvars(bsData{bin},Fly');
%     bsData{bin}.Properties.VariableNames{'Var4'} = 'Fly';
%     
%     mdl_bs_bm{bin} = fitlme(bsData{bin},'bump_var~1+bump_mag+(1|Fly)');
%     coef_bm(bin) = mdl_bs_bm{bin}.Coefficients.Estimate(2);
%     pval_bm(bin) = mdl_bs_bm{bin}.Coefficients.pValue(2);
%     
%     mdl_bs_hw{bin} = fitlme(bsData{bin},'bump_var~1+bump_width+(1|Fly)');
%     coef_hw(bin) = mdl_bs_hw{bin}.Coefficients.Estimate(2);
%     pval_hw(bin) = mdl_bs_hw{bin}.Coefficients.pValue(2);
%     
% end
% 
% %repeat figure with these coefficients
% figure('Position',[100,100,1400,800]),
% subplot(1,2,1)
% plot(bins,coef_bm,'color',[.5 .5 .5],'linewidth',2)
% hold on
% for bin = 1:length(bins)
%     if pval_bm(bin) < 0.05
%         plot(bins(bin),coef_bm(bin),'ro','linewidth',2)
%     else
%         plot(bins(bin),coef_bm(bin),'ko','linewidth',2)
%     end
%     hold on
% end
% ylim([-1 1]);
% title('Relationship between bump magnitude and bump instability (bump var)');
% ylabel('Regression coefficient');
% xlabel('Bin number');
% line([0 100],[0 0],'color','b','linewidth',2);
% 
% subplot(1,2,2)
% plot(bins,coef_hw,'color',[.5 .5 .5],'linewidth',2)
% hold on
% for bin = 1:length(bins)
%     if pval_hw(bin) < 0.05
%         plot(bins(bin),coef_hw(bin),'ro','linewidth',2)
%     else
%         plot(bins(bin),coef_hw(bin),'ko','linewidth',2)
%     end
%     hold on
% end
% ylim([-1 1]);
% title('Relationship between bump width and bump instability (bump var)');
% ylabel('Regression coefficient');
% xlabel('Bin number');
% line([0 100],[0 0],'color','b','linewidth',2);
% 
% %save figure
% saveas(gcf,[path,'\globalPlots\bumpStabilityThresh_acrossBlocks_pval_with_intercept.png']);
% 
% %% Look at the magnitude of the effect on bump variation
% 
% %look at the estimates in the model and how much they're affected by
% %changes in bump magnitude or bump width
% for bin = 1:length(bins)
%     bump_var_BM_zero(bin) = mdl_bs_bm{1,bin}.Coefficients.Estimate(1);
%     bump_var_BM_one(bin) = mdl_bs_bm{1,bin}.Coefficients.Estimate(1)+mdl_bs_bm{1,bin}.Coefficients.Estimate(2);
%     bump_var_BM_two(bin) = mdl_bs_bm{1,bin}.Coefficients.Estimate(1)+2*mdl_bs_bm{1,bin}.Coefficients.Estimate(2);
%     
%     bump_var_HW_two(bin) = mdl_bs_hw{1,bin}.Coefficients.Estimate(1)+2*mdl_bs_hw{1,bin}.Coefficients.Estimate(2);
%     bump_var_HW_three(bin) = mdl_bs_hw{1,bin}.Coefficients.Estimate(1)+3*mdl_bs_hw{1,bin}.Coefficients.Estimate(2);
% end
% 
% figure('Position',[100 100 1000 400]),
% %Plot the percentage change 
% subplot(1,2,1)
% change_BM = 100*(bump_var_BM_two-bump_var_BM_one)./bump_var_BM_one;
% plot(bins,change_BM,'-ko');
% hold on
% line([0 100],[0 0],'linewidth',2,'color','b');
% xlabel('Bin number');
% ylabel('Percentage change in bump var');
% title('Change in bump var per unit change in BM');
% 
% subplot(1,2,2)
% change_HW = 100*(bump_var_HW_three-bump_var_HW_two)./bump_var_HW_two;
% plot(bins,change_HW,'-ko');
% hold on
% line([0 100],[0 0],'linewidth',2,'color','b');
% title('Change in bump var per unit change in HW');
% xlabel('Bin number');
% 
% %save figure
% saveas(gcf,[path,'\globalPlots\bumpStability_acrossBlocks_effect_magnitude.png']);

%% Use zscored data

for bin = 1:length(bins)
    
    %Initialize empty table
    bsDataz{bin} = array2table(zeros(0,5),'VariableNames', {'bump_var','bump_mag','bump_width','zscoredBM','zscoredBW'});
    Fly = [];
    for fly = 1:length(data)
        Fly = [Fly,repelem(fly,length(data(fly).bs_table{1,bin}.bump_mag))];
        %Add zscored BM and BW per fly to the tables
        zscoredBM = zscore(data(fly).bs_table{1,bin}.bump_mag);
        data(fly).bs_table{1,bin} = addvars(data(fly).bs_table{1,bin},zscoredBM,'NewVariableNames','zscoredBM');
        zscoredBW = zscore(data(fly).bs_table{1,bin}.bump_width);
        data(fly).bs_table{1,bin} = addvars(data(fly).bs_table{1,bin},zscoredBW,'NewVariableNames','zscoredBW'); 
        %Combine all tables into one
        bsDataz{bin} = [bsDataz{bin};data(fly).bs_table{1,bin}];
    end
    %Add the fly ID as a variable
    bsDataz{bin} = addvars(bsDataz{bin},Fly');
    bsDataz{bin}.Properties.VariableNames{'Var6'} = 'Fly';
    
    %Fit a model to explain bump variability through bump mag or bump width
    mdl_bs_bm{bin} = fitlme(bsDataz{bin},'bump_var~1+bump_mag+(1|Fly)');
    coef_bm(bin) = mdl_bs_bm{bin}.Coefficients.Estimate(2);
    pval_bm(bin) = mdl_bs_bm{bin}.Coefficients.pValue(2);
    
    mdl_bs_hw{bin} = fitlme(bsDataz{bin},'bump_var~1+bump_width+(1|Fly)');
    coef_hw(bin) = mdl_bs_hw{bin}.Coefficients.Estimate(2);
    pval_hw(bin) = mdl_bs_hw{bin}.Coefficients.pValue(2);    
end

%% Plot the coefficients obtained from the model, with their significance

figure('Position',[100,100,1400,800]),
subplot(1,2,1)
plot(bins,coef_bm,'color',[.5 .5 .5],'linewidth',2)
hold on
for bin = 1:length(bins)
    if pval_bm(bin) < 0.05
        plot(bins(bin),coef_bm(bin),'ro','linewidth',2)
    else
        plot(bins(bin),coef_bm(bin),'ko','linewidth',2)
    end
    hold on
end
ylim([-1 1]);
title('Relationship between bump magnitude and bump instability (bump var)');
ylabel('Regression coefficient');
xlabel('Bin number');
line([0 100],[0 0],'color','b','linewidth',2);

subplot(1,2,2)
plot(bins,coef_hw,'color',[.5 .5 .5],'linewidth',2)
hold on
for bin = 1:length(bins)
    if pval_hw(bin) < 0.05
        plot(bins(bin),coef_hw(bin),'ro','linewidth',2)
    else
        plot(bins(bin),coef_hw(bin),'ko','linewidth',2)
    end
    hold on
end
ylim([-1 1]);
title('Relationship between bump width and bump instability (bump var)');
ylabel('Regression coefficient');
xlabel('Bin number');
line([0 100],[0 0],'color','b','linewidth',2);

%Save figure
saveas(gcf,[path,'\globalPlots\bumpStability_acrossBlocks.png']);

%% Compute and plot the magnitude of the effect

%We're going to quantify what percentage change in bump variability by a
%1 SD increase in bump magnitude or bump width (since the data is zscored,
%this accounts to comparing the bump variability estimate for BM or BW = 0
%vs BM or BM =1

for bin = 1:length(bins)
    bump_var_BM_zero(bin) = mdl_bs_bm{1,bin}.Coefficients.Estimate(1);
    bump_var_BM_one(bin) = mdl_bs_bm{1,bin}.Coefficients.Estimate(1)+mdl_bs_bm{1,bin}.Coefficients.Estimate(2);
    
    bump_var_HW_zero(bin) = mdl_bs_hw{1,bin}.Coefficients.Estimate(1);
    bump_var_HW_one(bin) = mdl_bs_hw{1,bin}.Coefficients.Estimate(1)+mdl_bs_hw{1,bin}.Coefficients.Estimate(2);
end

figure('Position',[100 100 1000 400]),
%Plot the percentage change 
subplot(1,2,1)
change_BM = 100*(bump_var_BM_one-bump_var_BM_zero)./bump_var_BM_zero;
plot(bins,change_BM,'-ko');
hold on
line([0 100],[0 0],'linewidth',2,'color','b');
xlabel('Bin number');
ylabel('Percentage change in bump var');
title('Change in bump var per 1 SD increase in BM');
ylim([-40 40]);

subplot(1,2,2)
change_HW = 100*(bump_var_HW_one-bump_var_HW_zero)./bump_var_HW_zero;
plot(bins,change_HW,'-ko');
hold on
line([0 100],[0 0],'linewidth',2,'color','b');
title('Change in bump var per 1 SD increase in HW');
xlabel('Bin number');
ylim([-40 40]);

%Save figure
saveas(gcf,[path,'\globalPlots\bumpStability_effectMagnitude.png']);

%% Import laser power control data for comparison

%list all of the files
lp_path = 'Z:\Wilson Lab\Mel\Experiments\Uncertainty\Laser_power_control\data';
folderNames = dir(lp_path);

for folder = 1:length(folderNames)
    if contains(folderNames(folder).name,'60D05')
        
        %load session info
        load(fullfile(folderNames(folder).folder,folderNames(folder).name,'sessions_info.mat'))
        
        %list the data file names
        fileNames = dir(fullfile(folderNames(folder).folder,folderNames(folder).name,'analysis'));
        
        %load the data corresponding to the different laser power values
        for file = 1:length(fileNames)
            if contains(fileNames(file).name,['sid_',num2str(sessions_info.laser_5)])
                lp_data{folder,1} = load(fullfile(fileNames(file).folder,fileNames(file).name));
            elseif contains(fileNames(file).name,['sid_',num2str(sessions_info.laser_10)])
                lp_data{folder,2} = load(fullfile(fileNames(file).folder,fileNames(file).name));
            elseif contains(fileNames(file).name,['sid_',num2str(sessions_info.laser_15)])
                lp_data{folder,3} = load(fullfile(fileNames(file).folder,fileNames(file).name));
            end
        end
    end  
end

%Remove empty cells
lp_data = lp_data(~cellfun(@isempty, lp_data(:,1)),:);

%% For the laser power control, compute the SD change in BM between conditions?

%Bin the data for each laser power
%Compute the mean zBM, zBW and variability per bin per laser power
%Compare the relationship between parameters across laser powers?

figure,

%bump magnitude
for fly = 1:length(lp_data)
    for laser_power = 1:3
        BM(fly,laser_power) = mean(lp_data{fly,laser_power}.data.bump_magnitude);
    end
end

subplot(1,3,1)
plot(BM','-o','color',[.5 .5 .5])
hold on
errorbar(1:3,mean(BM),std(BM)/sqrt(length(BM)),'-ko');
xlim([0 4]);
xticks([1 2 3]);
xticklabels({'5','10','15'});
xlabel('Laser power (%)');
ylabel('Bump magnitude');


%bump width
for fly = 1:length(lp_data)
    for laser_power = 1:3
        half_max_width = compute_bump_width(lp_data{fly,laser_power}.data.mean_dff_EB);
        BW(fly,laser_power) = mean(half_max_width);     
    end
end

subplot(1,3,2)
plot(BW','-o','color',[.5 .5 .5])
hold on
errorbar(1:3,mean(BW),std(BW)/sqrt(length(BW)),'-ko');
xlim([0 4]);
xticks([1 2 3]);
xticklabels({'5','10','15'});
xlabel('Laser power (%)');
ylabel('Bump width');


%bump_variability
for fly = 1:length(lp_data)
    for laser_power = 1:3
        bump_pos = lp_data{fly,laser_power}.data.dff_pva;
        [var bump_var(fly,laser_power)] = circ_std(bump_pos);
    end
end

subplot(1,3,3)
plot(bump_var','-o','color',[.5 .5 .5])
hold on
errorbar(1:3,mean(bump_var),std(bump_var)/sqrt(length(bump_var)),'-ko');
xlim([0 4]);
xticks([1 2 3]);
xticklabels({'5','10','15'});
xlabel('Laser power (%)');
ylabel('Bump variability');
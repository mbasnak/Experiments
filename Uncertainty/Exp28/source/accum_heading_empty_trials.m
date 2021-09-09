%Code to generate a figure comparing accumulated heading and bump tracking
%in empty trials

clear all; close all;
folderNames = dir('Z:\Wilson Lab\Mel\Experiments\Uncertainty\Exp28\data');


%for the fly folders
for folder = 1:length(folderNames)
    if (contains(folderNames(folder).name,'60D05')==1 & contains(folderNames(folder).name,'20210127')==0 & contains(folderNames(folder).name,'20210218_60D05_7f_fly2')==0) %I'm excluding the first fly because I didn't run a darkness trial for it
        
        path = [folderNames(folder).folder,'\',folderNames(folder).name];
        %get the sessions info
        load([path,'\sessions_info.mat'])
        
        %load the empty trial and inverted gain trial data
        load([path,'\analysis\analysis_sid_',num2str(sessions_info.empty_trial),'_tid_0.mat']);
               
        figure('Position',[200 200 1000 600]),
        fly_position = -unwrap(data.heading);
        plot(data.time,fly_position,'LineWidth',1.5,'color',[0.4660, 0.6740, 0.1880])
        hold on
        phase = unwrap(data.dff_pva);
        %shift phase to match fly_position
        shifted_phase = phase - abs(fly_position(1));
        plot(data.time,shifted_phase,'LineWidth',1.5,'color',[0.6350, 0.0780, 0.1840])
        xlim([0 120]);
        axP = get(gca,'Position');
        legend('Walking','Bump estimate','Location', 'Best')
        set(gca, 'Position', axP);
        ylabel('Accumulated heading (rad)','fontweight','bold','fontsize',12);
        xlabel('Time (sec)','fontweight','bold','fontsize',12);
        
        saveas(gcf,[path,'\analysis\plots\empty_trial_accum_heading.png']);     
        
    end
end









filePath = mfilename('fullpath')
folders = dir([filePath(1:end-26),'data\three_patterns']);

for i = 1:length(folders)
    if strfind(folders(i).name,'60D05') ~= 0
        files = dir([folders(i).folder,'\',folders(i).name]);
        
        for f = 1:length(files)
            if strfind(files(f).name,'analysis') ~= 0
        
                load([files(f).folder,'\',files(f).name])
        
                %Find frames corresponding to the three stim
                singleBar = find(data.fr_y_ds<4);
                ambiguous = find(data.fr_y_ds>6 & data.fr_y_ds<8);
                darkness = find(data.fr_y_ds>9);

                %with second heading and offset
                figure('Position',[100,100,1400,1000]),
                subplot(3,1,1)
                imagesc(data.dff_matrix)
                colormap(gray)
                ylabel('PB glomerulus');
                set(gca,'XTickLabel',[]);

                subplot(3,1,2)
                plot(data.time,rad2deg(data.flyPosRad))
                hold on
                plot(data.time,-rad2deg(data.phase))
                legend('fly heading', 'EPG phase');
                ylim([-180, 180]);
                xlim([0,data.time(end)]);
                testOffset = circ_dist(data.flyPosRad',data.phase);
                ylabel('PB glomerulus');
                set(gca,'XTickLabel',[]);

                %change color acording to stim
                subplot(3,1,3)
                plot(data.time(singleBar),rad2deg(testOffset(singleBar)),'.b')
                hold on
                plot(data.time(ambiguous),rad2deg(testOffset(ambiguous)),'.r')
                plot(data.time(darkness),rad2deg(testOffset(darkness)),'.k')
                xlim([0,data.time(end)]);
                ylim([-180, 180]);
                ylabel('Offset (deg)');
                xlabel('Time (sec)');
                legend('vertical bar', 'horizontal bar', 'panels off')
                
                saveas(gcf,[filePath(1:end-26),'plots\IndividualPlotFly',num2str(i),'.png']);
                close;

%% Get polar histograms for the heading and offset separated by stimulus tipe
                figure,
                subplot(3,2,1)
                polarhistogram(data.flyPosRad(singleBar),20,'FaceColor','b')
                title('Fly heading');
                subplot(3,2,2)
                polarhistogram(testOffset(singleBar),20,'FaceColor','b')
                title('Offset');

                subplot(3,2,3)
                polarhistogram(data.flyPosRad(ambiguous),20,'FaceColor','r')
                subplot(3,2,4)
                polarhistogram(testOffset(ambiguous),20,'FaceColor','r')

                subplot(3,2,5)
                polarhistogram(data.flyPosRad(darkness),20,'FaceColor','k')
                subplot(3,2,6)
                polarhistogram(testOffset(darkness),20,'FaceColor','k')
                
                saveas(gcf,[filePath(1:end-26),'plots\IndPolarHistFly',num2str(i),'.png']);
                close;
            end
        end

    end
end
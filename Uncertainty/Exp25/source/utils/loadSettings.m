function [NiDaqChannels] = loadSettings(room)

if strcmp(room,'Jenny-2P') == 1
    NiDaqChannels.headingFly = 1;
    NiDaqChannels.yawFly = 2;
    NiDaqChannels.xFly = 3;
    NiDaqChannels.xFlyGain = 4;
    NiDaqChannels.xPanels = 5;
    NiDaqChannels.yPanels = 6;
    NiDaqChannels.PanelStatus = 7;
    NiDaqChannels.OptoTrigger = 8;
        
elseif strcmp(room,'Jenny-Behavior') == 1    
    NiDaqChannels.headingFly = 1;
    NiDaqChannels.yFly = 2;
    NiDaqChannels.xFly = 3;
    NiDaqChannels.xFlyGain = 7;
    NiDaqChannels.xPanels = 4;
    NiDaqChannels.yPanels = 5;
    NiDaqChannels.PanelStatus = 6;
    NiDaqChannels.OptoTrigger = 8;
    
elseif strcmp(room,'Mel-2P') == 1
    NiDaqChannels.headingFly = 1;
    NiDaqChannels.yFly = 2;
    NiDaqChannels.xFly = 3;
    NiDaqChannels.xFlyGain = 7;
    NiDaqChannels.xPanels = 4; %x and y might be inverted - check (MB 20191005)
    NiDaqChannels.yPanels = 5;
    NiDaqChannels.PanelStatus = 6;
    NiDaqChannels.OptoTrigger = 8;   
    
else
    NiDaqChannels.headingFly = 1;
    %NiDaqChannels.yFly = 2; 
    NiDaqChannels.yawFly = 2; %changed on 20200110 to match new channel arrangement
    NiDaqChannels.xFly = 3;
    NiDaqChannels.xFlyGain = 7;
    NiDaqChannels.xPanels = 4;
    NiDaqChannels.yPanels = 5;
    NiDaqChannels.PanelStatus = 6;
    NiDaqChannels.OptoTrigger = 8;



end
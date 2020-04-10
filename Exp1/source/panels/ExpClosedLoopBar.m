% Code to run the NiDaq and a closed-loop bar and save the data

function ExpClosedLoopBar(flyNum,expNum,time,folder)

%INPUTS
%flyNum = what number of fly in the day is this one
%expNum = what exp number is it for this particular fly
%time = how long do you want to run the closed-loop bar for
%folder = what experiment number is this of all my experiments, to later
%save the data there.


% Move to the directory where I want to save the data to
cd(strcat('Z:\Wilson Lab\Mel\Experiments\Exp1',num2str(folder),'\data\'));

daqreset %reset DAC object
devID = 'Dev1';  % Set device ID (to know what the ID is, you can type "daq.getDevices"

% Configure session: national instruments output/input
niOI = daq.createSession('ni'); %create a session
niOI.Rate = 1000;% set sample rate
niOI.DurationInSeconds = time; %set duration in seconds
% Determine the analog INPUT Channels
aI = niOI.addAnalogInputChannel( devID , 1:7 , 'Voltage' );
% Set all channels to the correct inputType, likely 'SingleEnded'
for i = 1:7
    aI(i).InputType = 'SingleEnded';
end

pattern = 1;
Panel_com('set_pattern_id', pattern); %load the light stripe pattern
pause(0.01)
Panel_com('send_gain_bias', [10 0 0 0]);
pause(0.01)
Panel_com('set_mode', [3 0]); %set the x to be controlled by FicTrac and the Y to be open loop
Panel_com('start');

% Acquire data
rawData = niOI.startForeground(); %this will acquire the channels described above for the length of time also defined

pause(2);
% Stop the panels
Panel_com('stop');
Panel_com('all_off');
Panel_com('set_AO',[3 0]);


if flyNum == 1 && expNum == 1 %if it's the first fly and the first experiment
   mkdir ([date]) %make a folder with today's date
end


if expNum == 1 %if it's the first experiment for this fly
   cd(strcat('Z:\Wilson Lab\Mel\Experiments\Exp1',num2str(folder),'\data\',date)); %move to today's folder
   mkdir (strcat('flyNum',num2str(flyNum))) %inside that folder make a folder for this fly
   cd(strcat('Z:\Wilson Lab\Mel\Experiments\Exp1',num2str(folder),'\data\',date,'\flyNum',num2str(flyNum)));
   getFlyInfo() %get fly's details
else
   cd(strcat('Z:\Wilson Lab\Mel\Experiments\Exp1',num2str(folder),'\data\',date,'\flyNum',num2str(flyNum))); %otherwise move to this fly's folder
end


save(strcat('dataClosedLoopBar',num2str(expNum),'.mat'),'rawData','pattern'); %save the data 


end
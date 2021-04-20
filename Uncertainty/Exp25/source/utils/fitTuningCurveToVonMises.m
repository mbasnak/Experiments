function [vonMises, rescaledVonMises] = fitTuningCurveToVonMises( data , angles_rad )

%This is Yvette's function to fit data to a vonMises distribution

% FITTUNIGNCURVETOVONMISES fits a tuning curve to a von mises function
% 
% Used to analize open loop of closed loop tuning data to fit the curve
%
% This function requires the CircStat tool box: 
% https://www.mathworks.com/matlabcentral/fileexchange/10676-circular-statistics-toolbox-directional-statistics
%
%   INPUTS
%
%   data - response amplitude at various angle locations/ bar
%   positions
%
%   angles_rad - angles, in radians, where the amplitude measurements occured (must be an array with same
%   number of elements at 'data')
%
%
%   OUTPUT 
%   vonMises - fit von mises function
%   rescaledVonMises - fit von mises function rescaled to match data values
%   
%   Yvette Fisher 8/2018

% FUNCTION STEPS:

% 0) TODO make sure that data is in collumns
data = data(:);
angles_rad = angles_rad(:);

% check that data and angle_rad are same size
if( size( data ) ~= size( angles_rad ) )
    error(' ERROR: data and angle_rad must be the same length');
end

% check that angle_rad values are all positive
if( sum( angles_rad < 0 ) > 0 )
    error(' ERROR: angle_rad values must all be positive');
end


% 1) normalize incoming data (min subtract and max divide) - save min and
% max values as variable
dataMin = min( data );

% min subtract
dataMinSubtract  = data - dataMin;

dataMaxAfterMinSubtract = max( dataMinSubtract );

% max divide
normData = dataMinSubtract / dataMaxAfterMinSubtract;


% 2) find von mises function paramter that fit the normalized data

% find parameters of von mises that fit normData
[thetahat, kappa ] = circ_vmpar( angles_rad' , normData  );

% 3) make von mises function using found parameters
 [p, ~ ] = circ_vmpdf(angles_rad' , thetahat, kappa);

% 4) normalize new von mises function ( min subtract and max divide) 

pMinSubtract = p - min(p); % min subtract
p_norm = pMinSubtract / max( pMinSubtract ); % max divide

% 5) rescale von mises function to match orignial data
p_rescaled = ( p_norm * dataMaxAfterMinSubtract ) + dataMin;

% % 6) plot original data and the rescaledVonMises for user comparison
%  figure; set(gcf, 'Color', 'w'); 
%  plot( angles_rad, data); hold on;
%  plot( angles_rad, p_rescaled); hold on;
 %niceaxes; box off;

% 7) output values
vonMises = p;
rescaledVonMises = p_rescaled;

end


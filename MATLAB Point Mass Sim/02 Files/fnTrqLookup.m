function torque = fnTrqLookup(veh,speed,isRPM)

%  Creates an array of torques from given motor data.
% 
%  Inputs :
%       veh     : file name or structure containing vehicle parameters
%                 (e.g., torque, mass etc.).
%                 veh.Trq is a string containing the motor data reference 
%                 for the excel file which is the only parameter which will
%                 be used here
%       speed   : a 1D array of speeds to look up the torque
%       isRPM   : boolean value, whether speed is in rpm or m/s
%  
%  Outputs :
%       torque  : a 1D arry of torque values

if nargin < 3
    isRPM = 0;
end

% Initialise torque array
torque = nan(size(speed));

% veh.Trq is the name of the sheet to read from
if ischar(veh.Trq)
    xlfile = '03 Motor Data/Motor Data.xlsx';
    % First column should be rpm and second is torque
    xldata = readmatrix(xlfile,'sheet',veh.Trq,'ExpectedNumVariables',2);
    veh.Trq = rmmissing(xldata); % remove nan values
end

for i = 1:length(torque)
    if isscalar(veh.Trq) % If using constant torque value
        torque(i) = veh.Trq;
    else
        if ~isRPM
            % Convert speed in m/s to rpm
            speed_in_rpm = veh.SprocRatio * (speed(i)/veh.RollRad);
            speed_in_rpm = 60*speed_in_rpm/(2*pi);
        else
            speed_in_rpm = speed(i);
        end
        % Interpolate motor rpms using motor data
        torque(i) = interp1(veh.Trq(:,1),veh.Trq(:,2),speed_in_rpm,'linear',veh.Trq(end,2));
    end
end

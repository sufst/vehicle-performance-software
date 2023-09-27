function trqcurve = createTrqCurve(veh, BPlot)
 
%  Function looks up a torque vs rpm curve for the given motor.
%
%  Written by    : Hayley Brooks
%  Last modified : 21-08-2023
% 
%  Inputs :
%       veh     : file name or structure containing vehicle parameters
%                 (e.g., torque, mass etc.). 
%                 more specifically veh.Trq is a string with the name of
%                 the motor to look up the data for
%       BPlot : (optional) boolean to plot torque curve 
%  
%  Outputs :
%       trqcurve : a 2D array with 2 columns as rpm and torque (Nm)

if nargin < 2
    BPlot = 0;
end

% Create arrays for rpm and torque data
%veh.MaxRPM = max(veh.Trq(:, 1));
motor_rpm = (0:10:veh.MaxRPM);

torque = fnTrqLookup(veh, motor_rpm, 1);

%% Plot torque curve
if BPlot
    figure()
    plot(motor_rpm, torque)
    xlabel('RPM')
    ylabel('Torque (Nm)')
    title([veh.name ' Torque Curve'])
end

trqcurve = [reshape(motor_rpm,[],1), reshape(torque,[],1)];

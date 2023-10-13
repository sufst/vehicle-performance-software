function res = runLapSim(veh,track,options)
 
%  Overall function to run the lap sim. Calls other functions to create the
%  cars GGV and then 'solves' the lap using the given data.
%
%  Written by    : Will Doyle
%  Last modified : 21-08-2023
% 
%  Inputs :
%       veh     : file name or structure containing vehicle parameters
%                 (e.g., torque, mass etc.). if empty will use StagVII.mat.
%       track   : file name or structure containing track parameters (e.g., 
%                 distance and curvature). if empty will use
%                 FSUK_Sprint_MP.mat
%       options : (optional) structure containing various options to configure
%                 lapsim. See LapSimSolver.m for more details.
%  
%  Outputs :
%       res     : structure containing all the channels, metrics, input data 
%                 and log from the LapSim
%
%  Example : res = runLapSim('StagVIII','FSUK_Sprint_MP')

%% Deal with inputs
addpath(genpath(pwd)); % add subfolders to path

% Import vehicle data
if nargin < 1 || isempty(veh)
    veh = load('01 Vehicles/StagX.mat');
    veh.Trq = 'Emrax 228 LC';
    warning('No vehicle specified, using StagX.mat')
elseif ischar(veh)
    veh = load(veh);
end

veh.GFLat = 0.6;
veh.GFLong = 0.6;

% Import track data
% Will not work with skip pad as no minimum curvature is found
% For accel want to use peak torque
if nargin < 2 || isempty(track)
    track = load('00 Tracks/FSUK_Sprint_MP.mat');
    warning('No track specified, using FSUK_Sprint_MP.mat')
elseif ischar(track)
    track = load(track);
end

if nargin < 3 || ~isfield(options, 'PlotGGV')
    options.PlotGGV = 0;
end

if nargin < 3 || ~isfield(options, 'PlotTrqCurve')
    options.PlotTrqCurve = 0;
end

%% Call main functions
veh.Trq = createTrqCurve(veh,options.PlotTrqCurve);
GGV = CreatePointMassGGV(veh,options.PlotGGV);
res = LapSimSolver(GGV,track,veh,options);

%% Uncomment to calculate score for a certain event
score = calcSprintScore(res.metrics.LapTime);
disp(score)
%score = calcAccelScore(res.metrics.LapTime);

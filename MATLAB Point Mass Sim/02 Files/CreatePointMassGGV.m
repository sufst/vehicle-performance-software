function GGV = CreatePointMassGGV(veh, BPlot)

%  Function creates a GGV envelope for the given vehicle. GGV envelope
%  consists of the maximum lateral and longitudinal forces at any given
%  speed.
%
%  Written by    : Will Doyle
%  Last modified : 21-08-2023
%
%  Inputs :
%       veh     : file name or structure containing vehicle parameters
%                 (e.g., torque, mass etc.)
%       BPlot   : (optional) boolean to plot GGV in figure
%  
%  Outputs :
%       GGV     : structure containing details of the GGV envelope

%% Set up
if nargin < 2
    BPlot = 1;
end

% Define global parameters
rho = 1.225;
g = 9.81;
vSteps = 30; % Number of speeds to calc acceleration at
CombinedSteps = 15; % Number of gLat points to calc per speed (2 gLongs per gLat)

%% Load in Vehicle Data
% ------------- Torque/Motor Modelling -------------
% If torque/rpm curve hasn't been loaded already do so here
if ischar(veh.Trq)
    veh.Trq = createTrqCurve(veh);
end

% Find maximum force from max torque
maxforce = veh.SprocRatio * max(veh.Trq(:,2)) / veh.RollRad;

% Find maximum speed of the vehicle
vMaxGearing = veh.RollRad * (veh.MaxRPM * (2*pi/60) / veh.SprocRatio);
vMaxMotor = sqrt((2*maxforce) / (rho*veh.CdA));
vMax = min(vMaxGearing,vMaxMotor);

% ------------- Tyre Data -------------
% Load in tyre coefficients
veh.MuY = loadtyredata(['04 Tyres/Lookup_tables/NFY/' veh.Tyres{1}]);
veh.MuX = loadtyredata(['04 Tyres/Lookup_tables/NFX/' veh.Tyres{2}]);
% In this point mass model, 1 tyre supports the load of the whole car so scale up
veh.MuY(1,:) = veh.MuY(1,:)*4;
veh.MuX(1,:) = veh.MuX(1,:)*4;

% Multiply tyre coefficient by grip factors
veh.MuY(2,:) = veh.MuY(2,:) * veh.GFLat;
veh.MuX(2,:) = veh.MuX(2,:) * veh.GFLong;

%% Initialise arrays
Speeds = [0 1 2 3 logspace(log10(4),log10(vMax),vSteps-4)]; % Speeds at which we will evaluate GG envelope. Logspace gives better spacing for curvature interpolation
gLat = zeros(vSteps,CombinedSteps);
gLongNeg = zeros(vSteps,CombinedSteps);
gLongPos = zeros(vSteps,CombinedSteps);

%% Calculate gLat and gLong for speed
for i = 1:length(Speeds)

    v = Speeds(i);
    
    % Calculate body forces on car
    BodyForces(i).Speed = v;
    BodyForces(i).FDrag = 0.5 * rho * veh.CdA * (v^2);
    BodyForces(i).FLift = 0.5 * rho * veh.ClA * (v^2);
    BodyForces(i).FGravity = veh.mass * g;
    % Calculate vertical force on tyre
    BodyForces(i).Fz = BodyForces(i).FLift + BodyForces(i).FGravity;
    
    %% GG Cases
    % Calc max lat accel case
    gLat(i,end) = CalcLatAcc(veh,BodyForces(i),0);
    
    % Loop over all possible lateral cases calculating max possible long
    % accels (+ve and -ve)
    gLat(i,:) = linspace(0,gLat(i,end),CombinedSteps);
    
    for j = 1:length(gLat(i,:))-1
        [gLongNeg(i,j), gLongPos(i,j)] = CalcLongAcc(veh, BodyForces(i), gLat(i,j));
    end
end


%% Process in to GGV table
GGV.Speed = repmat(Speeds', 1, (2*CombinedSteps) - 1);
GGV.gLat = [gLat fliplr(gLat(:,1:end-1))];
GGV.gLong = [gLongNeg, fliplr(gLongPos(:,1:end-1))];

% Include other useful channels
GGV.Curvature = (GGV.gLat * g) ./ (GGV.Speed.^2);
GGV.Curvature(isnan(GGV.Curvature)) = 0;
GGV.Radius = 1./GGV.Curvature;

% Plot GGV diagram
if BPlot
    figure()
    mesh(GGV.gLong,GGV.gLat,GGV.Speed)
    xlabel('gLong')
    ylabel('gLat')
    zlabel('Speed (m/s)')
    title([veh.name ' GGV Envelope'])
end
end


%% %%%%%%%%%%%%%%%%%%%%%%% Functions %%%%%%%%%%%%%%%%%%%%%%% %%

function gLat = CalcLatAcc(veh, BodyForces, gLong)

g = 9.81;

% Calculate pure forces available from tyre
muY = interp1(veh.MuY(1,:), veh.MuY(2,:), BodyForces.Fz,'linear','extrap');
FyMax = muY * BodyForces.Fz;

muX = interp1(veh.MuX(1,:), veh.MuX(2,:), BodyForces.Fz,'linear','extrap');
FxMax = muX * BodyForces.Fz;

% Find the split between long/lat
Fx = (gLong * g * veh.mass) + BodyForces.FDrag;
Fy = solveFrictionShape(FxMax, FyMax, Fx, veh.GGVFactor);

gLat = Fy / (veh.mass * g);

end


function [gLongNeg, gLongPos] = CalcLongAcc(veh, BodyForces, gLat)

g = 9.81;

% Calculate pure forces available from tyre
muY = interp1(veh.MuY(1,:), veh.MuY(2,:), BodyForces.Fz,'linear','extrap');
FyMax = muY * BodyForces.Fz;

muX = interp1(veh.MuX(1,:), veh.MuX(2,:), BodyForces.Fz,'linear','extrap');
FxMax = muX * BodyForces.Fz;
% Find the split between long/lat
Fy = (gLat * g * veh.mass);
Fx = solveFrictionShape(FyMax, FxMax, Fy, veh.GGVFactor);

% Calculate how much force is available from PU
Trq = fnTrqLookup(veh,BodyForces.Speed);
FPU = veh.SprocRatio * veh.DriveEfficiency * Trq / veh.RollRad;

% Add drag (+ Rolling Resistance at some point?)
FLongNeg = -Fx - BodyForces.FDrag;
FLongPos = min(FPU, veh.WDist*Fx) - BodyForces.FDrag; % WDist is factor to represent only rear wheels can provide +ve long force

gLongNeg = FLongNeg / (veh.mass * g);
gLongPos = FLongPos / (veh.mass * g);

end


function F2 = solveFrictionShape(F1Max, F2Max, F1, shapefactor)
% Takes "pure" forces in each direction and tyre forced used in 1 direction to
% return tyre force available in other direction. 

%         FCirc = sqrt(Ftotal^2 - F1^2);
        FElipse = F2Max * sqrt(1 - (F1^2/F1Max^2)); 
%         FDiam = Ftotal - F1;
        FDiam = F2Max - (F2Max/F1Max)*F1;
        F2 = (shapefactor * FElipse) + ((1-shapefactor) * FDiam);
end


function data = loadtyredata(file)
d = load(file);
flds = fields(d);
if length(flds) ~= 1
    error('Cannot parse tyre data from %s\nmust include only 1 variable',file)
end
data = d.(flds{1});
end
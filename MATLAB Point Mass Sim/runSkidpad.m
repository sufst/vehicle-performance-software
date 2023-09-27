function res = runSkidpad(veh,r)

%  Function solves a lap sim for a skidpad event based on a given car and
%  radius of the circles.
%
%  Inputs:
%       veh     : structure containing vehicle data
%       r       : radius - default = 9.125m
%
%  Outputs :
%       res     : structure containing all the metrics from the simulation


if nargin < 1 || isempty(veh)
    veh = load('01 Vehicles/StagVIII.mat');
    warning('No vehicle specified, using StagVIII.mat')
end

if nargin < 2
    r = 9.125;
end

addpath('02 Files')
curv = 1/r;
dist = 2*pi*r;

ggv = CreatePointMassGGV(veh,0);

iPureLat = find(ggv.gLong(1,:)==0);
BGood = isfinite(ggv.Curvature(:,iPureLat));

res.speed = interp1(ggv.Curvature(BGood,iPureLat), ggv.Speed(BGood,iPureLat), curv);

res.time = dist/res.speed;
res.vkmh = res.speed * 3.6;
res.dist = dist;
res.gLat = ((res.speed^2) / r) / 9.81;

% Score will typically be 75 as this is ideal.
%res.score = calcSkidScore(res.time);

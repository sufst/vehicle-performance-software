function [x,y] = fnCoordPlot(track,BPlot)
% Input name of track definition file, output x,y coords of track map
% Author : Kris Collier, Marwan Ismail, Will Doyle 11/20

if nargin < 2
   BPlot = 1; 
end

if ischar(track)
    track = load(track);
end

Radius = 1./track.Curvature;

x = zeros(size(track.Curvature));
y = x;
heading = 0;

for i=1:length(track.Curvature)
    if track.Curvature(i) % if it's non-zero
        theta = track.Curvature(i);
        radius = Radius(i);
        
        thetaDivTwo = theta / 2;
        chord = 2*radius*sin(thetaDivTwo);

    else % maths falls apart for 0 curvature, hard code values
        theta = 0;
        chord = 1;
        thetaDivTwo = 0;
    end

    cCosTheta = chord*(cos(heading + thetaDivTwo));
    cSinTheta = chord*(sin(heading + thetaDivTwo));
    
    heading = heading + theta;
    
    x(i+1) = x(i) + cCosTheta;
    y(i+1) = y(i) + cSinTheta;

end
if BPlot
    plot(x,y,'LineWidth',5)
    axis equal
end
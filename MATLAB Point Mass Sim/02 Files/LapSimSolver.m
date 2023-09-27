function out = LapSimSolver(GGV, track, veh, opts)

%  Takes a Vehicle's GGV envelope and a track trajectory, and solves the 
%  fastest possible lap, operating on surface of GGV.
%  A quasi static lap sim solver using a forwards/reverse method.
% test
%  Written by    : Will Doyle
%  Last modified : 22-08-2023
% 
%  Inputs :
%       GGV     : structure containing data points of speed, gLat, gLong,
%                 curvature that defines the vehicles performance envelope
%       track   : structure containing curvature (1/m) + distance data at 
%                 1 metre intervals
%       veh     : vehicle structure defining car parameters
%                 used to post process channels
%       opts    : (optional) structure containing any options you may want to
%                 change from default

%                 StartSpeed     : default = 0
%                                  Can 'hard code' speed at start line, (0 for standing start) or set to -1 for cyclic speed
%                 BPlotProgress  : default = 0 
%                                  Displays live updates of speed calcs, useful for debuging but slows the program down
%                 BOutputLog     : default = 1
%                                  Outputs log messages to command window
%                 tol            : default = 10%  
%                                  Percentage speed by which solver is allowed to go past previously calc'd theoretical max speed. Aids in slow corners.
%                 LiftOffPercent : default = 0
%                                  Percentage of end of straight liftoff. Explained further in postproLiftOff.m.
%  
%  Outputs :
%       out     : structure containing all the channels, metrics, input data 
%                 and log from the LapSim


%% Set up
tic % Start timer
% Check if any options were input otherwise create default
optdefaults = { 'StartSpeed',       0;
                'BPlotProgress',    0;
                'BOutputLog',       1;
                'tol',              10;
                'LiftOffPercent',   0;
                };
if nargin < 3
    opts = struct();
end

for i = 1:length(optdefaults)
    if ~isfield(opts, optdefaults{i,1})
        opts.(optdefaults{i,1}) = optdefaults{i,2};
    end
end    

% Create log object to contain messages output to cmd window
messageLog = OutputLogger(opts.BOutputLog);


%% Find max possible speed at each point around the track then find apices
messageLog.addMsg('Calculating max possible speed at every point')

% Get maximum possible vehicle speed from GGV diagram
vMax = max(max(GGV.Speed));

% Create array of maximum speeds the length of the track data
MaxSpeed = vMax * ones(1,length(track.Curvature));

% Find the minimum curvature -> closer to straight line
% Should be min(min(?
min_curv = max(min(GGV.Curvature));

for i = 1:length([track.Curvature])
    % If curvature is bigger than minimum, find max speed the curve can be
    % taken, else use max vehicle speed
    if abs(track.Curvature(i)) > min_curv
        % Linearly interpolate maximum speeds from GGV data
        MaxSpeed(i) = fn2DInterp(GGV.gLong,GGV.Curvature,GGV.Speed,0,abs(track.Curvature(i)));
    end
end

% Find local minimum points i.e., apecies
% BApex is an array of 1s and 0s -> 1s indicate an apex
BApex = islocalmin(MaxSpeed);
% Array of indexes where apecies are
idApex = find(BApex);
% Create apex structure
Apex.Distance = [];
Apex.Curvature = [];
Apex.Speed = [];
for i = 1:length(idApex)
    Apex(i).Distance = track.Distance(idApex(i));
    Apex(i).Curvature = track.Curvature(idApex(i));
    % Find apex speed from GGV data
    Apex(i).Speed = fn2DInterp(GGV.gLong,GGV.Curvature,GGV.Speed,0,abs(Apex(i).Curvature));
    messageLog.addMsg(sprintf('Found apex at %4.1f m. Max Speed = %4.1f m/s',Apex(i).Distance,Apex(i).Speed));
end

[~,I] = sort([Apex.Speed]);
Apex = Apex(I);
messageLog.addMsg(sprintf('%d Apices found in total',i));

if opts.BPlotProgress
    figure()
    s(1) = subplot(2,1,1); hold on;
    pltLines.MaxSpeed = plot(track.Distance,MaxSpeed,'g');
    pltLines.Apices = plot([Apex.Distance],[Apex.Speed],'g.','MarkerSize',12);
    ylabel('Speed (m/s)')
    
    s(2) = subplot(2,1,2); hold on;
    plot(track.Distance,track.Curvature,'k')
    plot([Apex.Distance],[Apex.Curvature],'g.','MarkerSize',12)
    ylabel('Curvature (1/m)')
    
    linkaxes(s,'x')
end

%% Add "Apex" at start
Apex(end+1).Distance = 0;
Apex(end).Curvature = 0;
Apex(end).Speed = opts.StartSpeed;
% if standing start, curvature MUST be zero for GGV interpolation
if opts.StartSpeed == 0
    track.Curvature(1:2) = 0;
end

%% Set up braking and accel zone calcs
% Split GGV so that only +ve/-ve half is used for lookup
[~,iNeutral] = min(abs(GGV.gLong(1,:)));
GFields = fields(GGV);
for i = 1:length(GFields)
    BrkGGV.(GFields{i}) = GGV.(GFields{i})(:,1:iNeutral);
    AccGGV.(GFields{i}) = GGV.(GFields{i})(:,iNeutral:end);
end

% Initialise arrays
BrakingSpeed = nan(size(MaxSpeed));
AccelSpeed = BrakingSpeed;
longAcc = zeros(size(MaxSpeed));

%% Braking zone
for i = 1:length(Apex)
    currdist = Apex(i).Distance;
    currspeed = Apex(i).Speed;
    
    BOK = 1;
    messageLog.addMsg(sprintf('Calculating braking zone for apex at %d m', currdist))
    
    while BOK && currdist > 0
        currID = currdist+1; % Matlab indexing starts at 1 not 0 haha funny meme
        curv = abs(track.Curvature(currID-1)); % Curvature of preceeding section
        % Calc what long acceleration (and thus speed) is possible given
        % current speed and curvature
        [newspeed,longAcc(currID)] = fnfindspeed(BrkGGV,currspeed,curv);
        
        % Check we've successfully calc'd a speed
        if isnan(newspeed)
            BOK = 0;
            messageLog.addMsg(sprintf('Solver returned NaN speed at %d m',currdist));
        end
        
        % Check we haven't already calc'd a slower max speed from previous
        % apex
        if ~isnan(BrakingSpeed(currID-1)) && newspeed > BrakingSpeed(currID-1)
            messageLog.addMsg(sprintf('Skipping to next apex as speed at %d m (%4.1f) is faster than existing speed (%4.1f)',track.Distance(currID-1), currspeed, BrakingSpeed(currID-1)));
            BOK = 0;
            continue
        end
        % Check we aren't now quicker than max speed possible
        if newspeed > MaxSpeed(currID-1) * (1+(opts.tol/100)) || newspeed > vMax
            messageLog.addMsg(sprintf('Intersected Max speed line at %d m (%.1f m/s vs %.1f m/s)',track.Distance(currID-1),newspeed,MaxSpeed(currID-1)))
            BOK = 0;
            continue
        end
        
        % Otherwise add to array and shift back to preceeding section
        BrakingSpeed(currID-1) = newspeed;
        currspeed = newspeed;
        currdist = currdist - 1;
        
        % Plot if applicable
        if opts.BPlotProgress
            subplot(2,1,1)
            if isfield(pltLines,'Braking')
                pltLines.Braking.YData = BrakingSpeed;
            else
                pltLines.Braking = plot(track.Distance,BrakingSpeed,'r.');
            end
            drawnow()
        end
    end
end

longAcc = longAcc * -1; % all the accelerations came out +ve

%% Acceleration zone
for i = 1:length(Apex)
    currdist = Apex(i).Distance;
    currspeed = Apex(i).Speed;
    
    BOK = 1;
    messageLog.addMsg(sprintf('Calculating acceleration zone for apex at %d m', currdist))
    
    while BOK
        currID = currdist+1;
        curv = abs(track.Curvature(currID+1)); % Curvature of upcoming section
        
        [newspeed,longAcc(currID)] = fnfindspeed(AccGGV,currspeed,curv);

        % Check we've successfully calc'd a speed
        if isnan(newspeed)
            BOK = 0;
            messageLog.addMsg(sprintf('Solver returned NaN speed at %d m',currdist));
        end
        % Check we haven't already calc'd a slower max speed accelerating 
        % from previous apex
        if ~isnan(AccelSpeed(currID+1)) && newspeed > AccelSpeed(currID+1)
            messageLog.addMsg(sprintf('Skipping to next apex as speed at %d m (%4.1f) is faster than existing speed (%4.1f)',track.Distance(currID+1),currspeed,AccelSpeed(currID+1)));
            BOK = 0;
            continue
        end
        % Check we aren't now quicker than max speed possible. If we are,
        % adjust our projected speed.
        if newspeed > MaxSpeed(currID+1)
            newspeed = MaxSpeed(currID+1);
            longAcc(currID) = (MaxSpeed(currID+1)^2-MaxSpeed(currID)^2)/(2*9.81);
        end

        
        % Finally, check we haven't intersected a braking zone
        if ~isnan(BrakingSpeed(currID+1)) && BrakingSpeed(currID+1) < newspeed
            messageLog.addMsg(sprintf('Intersected Braking Zone at %d m',currdist))
            BOK = 0;
            continue
        end
        
        
        % Otherwise add to array and shift on to next section
        AccelSpeed(currID+1) = newspeed;
        currspeed = newspeed;
        currdist = currdist + 1;
        
        % Check if we've reached EOL. Time to sort of Start Of Lap if so
        if currdist == max(track.Distance)
            messageLog.addMsg('End of Lap reached')
            if opts.StartSpeed == -1
                Apex(end).Speed = currspeed; % cyclic lap speeds
                MaxSpeed(1) = currspeed;
            else
                MaxSpeed(1) = opts.StartSpeed; % Set this for the postpro later
            end
            BOK = 0;
        end
        
        % Plot if applicable
        if opts.BPlotProgress
            subplot(2,1,1)
            if isfield(pltLines,'Accel')
                pltLines.Accel.YData = AccelSpeed;
            else
                pltLines.Accel = plot(track.Distance,AccelSpeed,'b.');
            end
            drawnow()
        end
    end
end

%% Post Processing
speedtrace = min([reshape(MaxSpeed,1,[]);reshape(BrakingSpeed,1,[]);reshape(AccelSpeed,1,[])]);
dist = reshape(track.Distance,1,[]);
% -------------------- Lift Off 
if opts.LiftOffPercent
    messageLog.addMsg('------ Calculating lift Off ------')
    resnoLO.speed = speedtrace;
    resnoLO.dist = dist;
    resLO = postproLiftOff(resnoLO,veh,opts.LiftOffPercent,messageLog);
    speedtrace = resLO.speed;
    dist = resLO.dist;
end

messageLog.addMsg('========== Complete ===========')

% ------------------------- Channels -------------------------
% ----------------------- Basic
out.channels.speed =  speedtrace;
out.channels.vkmh = out.channels.speed * 3.6;
out.channels.gLongRAW = reshape(longAcc,1,[]);
out.channels.dist = dist;
out.channels.curvature = reshape(track.Curvature,1,[]);

% ----------------------- Calculated
out.channels.time = cumsum([0 out.channels.speed(2:end).^-1]);

out.channels.gLat = (out.channels.curvature .* (out.channels.speed.^2))/9.81;
out.channels.gLong = [diff(out.channels.speed) ./ diff(out.channels.time) 0]/9.81; 

out.channels.FDrag = 0.5 * 1.225 * veh.CdA * (out.channels.speed.^2);
out.channels.FLift = 0.5 * 1.225 * veh.ClA * (out.channels.speed.^2);

out.channels.FTractive = ((9.81 * out.channels.gLong) * veh.mass) + out.channels.FDrag;
out.channels.FTractPos = out.channels.FTractive .* (out.channels.FTractive>0);
out.channels.FTractNeg = -1 * out.channels.FTractive .* (out.channels.FTractive<0);
out.channels.PTractive = out.channels.FTractive .* out.channels.speed;
out.channels.PTractPos = out.channels.FTractPos .* out.channels.speed;
out.channels.PTractNeg = out.channels.FTractNeg .* out.channels.speed;

out.channels.ETractive = cumtrapz(out.channels.time, out.channels.PTractive) / 1e6;
out.channels.ETractPos = cumtrapz(out.channels.time, out.channels.PTractPos) / 1e6;
out.channels.ETractNeg = cumtrapz(out.channels.time, out.channels.PTractNeg) / 1e6;

out.channels.TrqPUAvail = fnTrqLookup(veh,out.channels.speed);
out.channels.TrqPUUsed = (out.channels.FTractPos * veh.RollRad)/(veh.SprocRatio*veh.DriveEfficiency);

out.channels.RPM = (30/pi)*(veh.SprocRatio .* out.channels.speed / veh.RollRad);
if length(out.channels.RPM) ~= length(out.channels.dist) % if end lap on shift then wont auto fill last bit of array
    out.channels.RPM(length(out.channels.dist)) = 0;
end
out.channels.PPU = out.channels.TrqPUUsed .* out.channels.RPM * (pi/30); % This will only be right if gearing isn't used

out.channels.EPU = cumtrapz(out.channels.time, out.channels.PPU) / 1e6;

% out.channels.Throttle = out.channels.PTractPos / veh.PPU;
out.channels.Throttle = out.channels.TrqPUUsed ./ out.channels.TrqPUAvail;
out.channels.Throttle(out.channels.Throttle>1) = 1;

% ------------------------- Metrics ------------------------- 
out.metrics.LapTime = max(out.channels.time);
out.metrics.SolutionTime = toc;
out.metrics.PUEnergyOUT_MJ = max(out.channels.EPU);
% out.metrics.NetEnergy_MJ = out.channels.ETractive(end);
out.metrics.AvgThrottle = mean(out.channels.Throttle);
out.metrics.vMax = max(out.channels.speed);

% Print Metrics
mfields = fields(out.metrics);
for i = 1:length(mfields)
    messageLog.addMsg(sprintf('%s : %.2f',mfields{i},out.metrics.(mfields{i})))
end

messageLog.trim();
out.log = messageLog.Messages;

out.inputs.GGV = GGV;
out.inputs.options = opts;
out.inputs.track = track;
out.inputs.veh = veh;
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Functions
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function vq = fn2DInterp(X,Y,V,xq,yq)
% Matlab doesnt have built in function to interpolate data in form we
% have :(
BGood = isfinite(X) & isfinite(Y) & isfinite(V);
x = reshape(X(BGood),[],1);
y = reshape(Y(BGood),[],1);
v = reshape(V(BGood),[],1);

vq = griddata(x,y,v,xq,yq);

% if isnan(vq)
%     if xq > max(x) || yq > max(y)
%         vq = 0;
%     end
% end
end
%%
function [v_out,longAcc] = fnfindspeed(GGV, v_in, curv)
% Interpolate GGV to find long acc possible for given speed+curvature
if isnan(v_in)
    keyboard
end
GGV.gLong = abs(GGV.gLong); % braking acceleration is postive because we go backwards in time
longAcc = fn2DInterp(GGV.Speed,GGV.Curvature,GGV.gLong,v_in,curv);
% SUVAT to find speed (v^2 = u^2 + 2as)
v_out = sqrt((v_in^2) + (2*9.81*longAcc));
end
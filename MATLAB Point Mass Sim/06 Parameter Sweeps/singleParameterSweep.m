function data = singleParameterSweep(veh, track, param, minimum, maximum, interval)

% create structure to output data
data = struct;
options.BOutputLog = 0; % veh and track needs to be loaded before

% go through every interval between the definied min and max
for i = minimum:interval:maximum
    % set parameter and run lap sim
    veh.(param) = i;
    result = runLapSim(veh,track, options);
    % get lap time from result of the lapsim
    time = result.metrics.LapTime;

    % add values to data variable
    if i == minimum
        data.(param) = i;
        data.Time = time;
    else
        data.(param)(end+1) = i;
        data.Time(end+1) = time;
    end
end

% plot acceleration time vs input parameter
plt = plot(data.(param), data.Time, '-o');
xlabel(param)
ylabel('Time (s)')
set(plt, 'Color', 'Black', 'LineWidth', 1.5)

% calculate minimum lap time and corresponding parameter value
[t, index] = min(data.Time);
val = data.(param)(index);

% display optimal values
disp("Minimum Time: " + t)
disp(param + ": " + val)

end
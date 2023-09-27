function score = calcEfficiencyScore(time, energy, laps)
%{
function calculates the FSUK efficiency score for an electric car

efficiency:
    for an electric car the energy conversion is 0.45kgCO2 per kWh
    teams are only considered if:
        1. fuel consumption doesn't exceed the equivalent of 60.06 kgCO2/100km
        2. lap time is under 145% of the fastest team that completed the event

inputs:
    time    : total endurance lap time
    energy  : energy consumed across endurance in kWh
    laps    : no. laps completed
%}

% numbers taken from FSUK 2023
T_min_per_lap = 58.306; % minimum averaged endurance time
CO2_min_per_lap = 0.0328; % minimum mass of CO2 per lap
effFactor_min = 0.046; % minimum possible efficiency factor
effFactor_max = 0.699; % maximum possible efficiency


% calculate time and CO2 per lap of the team being scored
CO2_per_lap = (energy * 0.45) / laps;
T_per_lap = time / laps;


% efficiency factor and score using the formulas given in FSUK 2023 rules
effFactor = (T_min_per_lap / T_per_lap) * (CO2_min_per_lap / CO2_per_lap);

score = 100 * (((effFactor_min / effFactor) - 1) / ...
    ((effFactor_min / effFactor_max) - 1));

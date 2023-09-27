function score = calcAccelScore(time, besttime)
%{
function calculates the FSUK score for the acceleration event from the lap 
time given

accel:
    maximum of 75 points
    5 points are awarded for finishing at least 1 run
    5 points are awarded for placing first in the top six runoff

inputs:
  time        : accel time + penalties
  besttime    : fastest vehicle time + penalties
%}

% if besttime isn't given, use the value from FSUK 2023
if nargin < 2
    besttime = 4.281;
end

%% Catch edge cases
if time <= besttime
    score = 75;
    warning("Fastest vehicle time, maximum points are awarded.")
    return
end

if time > (besttime * 1.50)
    score = 5;
    warning("Time is outside the maximum time (150% of best time). 5 points are awarded for completing a run.")
    return
end

%% Calculate the score using the formula given in the FSUK 2023 rules
score = (65 * ((((besttime * 1.50) / time) - 1) / ...
   (((besttime * 1.50) / besttime) - 1))) + 5;

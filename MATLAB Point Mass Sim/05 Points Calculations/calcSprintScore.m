function score = calcSprintScore(time, besttime)
%{
function calculates the FSUK score for the sprint from the lap time given

sprint:
    maximum of 100 points
    5 points are awarded for finishing at least 1 run

inputs:
  time        : lap time to get the score for + penalties
  besttime    : the best lap time of the competition + penalties
%}

% if besttime isn't given, use the value from FSUK 2023
if nargin < 2
    besttime = 63.692;
end

%% Catch edge cases
if time <= besttime
    score = 100;
    warning("Fastest vehicle time, maximum points are awarded.")
    return
end

if time > (besttime * 1.50)
    score = 5;
    warning("Time is outside the maximum time (150% of best time). 5 points are awarded for completing a run.")
    return
end

% calculate the score using the formula given in the FSUK 2023 rules
score = (95 * ((((besttime * 1.45) / time) - 1) / ...
        (((besttime * 1.45) / besttime) - 1))) + 5;
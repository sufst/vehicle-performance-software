function score = calcSkidScore(time, besttime)
%{
function calculates the FSUK score for skidpad from the lap time given

skidpad:
    maximum of 75 points
    5 points are awarded for finishing at least 1 run

inputs:
  time        : average of timed left and right circle + penalties
  besttime    : the best lap time of the competition + penalties
%}

% if besttime isn't given, use the value from FSUK 2023
if nargin < 2
    besttime = 5.806;
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

% calculate the score using the formula given in the FSUK 2023 rules
tRatioTop = ((besttime * 1.25) / time)^2 - 1;
tRatioBottom = ((besttime * 1.25) / besttime)^2 - 1;
score = (70 * (tRatioTop / tRatioBottom)) + 5;
function score = calcEnduroScore(time, besttime)
%{
function calculates the FSUK score for endurance from the lap time given

endurance:
    maximum of 250 points
    25 points are awarded for finishing

inputs:
  time        : endurance time + penalties
  besttime    : fastest vehicle time + penalties
%}

% if besttime isn't given, use the value from FSUK 2023
if nargin < 2
    besttime = 1586.250;
end

%% Catch edge cases
if time <= besttime
    score = 250;
    warning("Fastest vehicle time, maximum points are awarded.")
    return
end

if time > (besttime * 1.50)
    score = 25;
    warning("Time is outside the maximum time (150% of best time). 5 points are awarded for completing a run.")
    return
end

% calculate the score using the formula given in the FSUK 2023 rules
score = (225 * ((((besttime * 1.45) / time) - 1) / ... 
    (((besttime * 1.45) / besttime) - 1))) + 25;
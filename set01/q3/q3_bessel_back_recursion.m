% 2/10/15 Stefan Countryman

% A script for calculating bessel functions using back recursion
% Based heavily on scripts by RDM

close all
clear all

% Pick arbitrary initial values for the recursion, since the
% unwanted eigenvalue will decay under back-recursion anyway

ultimateTerm = 1;
penultimateTerm = 0;

% Specify array of desired x values

X = 1:10;

% Specify which Bessel function to calculate (must be at least 1 because
% MATLAB inexplicably likes to start counting at 1)

N = 2;

% Specify an array of desired decimal precision values

P = [4 8];

% Specify target accuracy

A = zeros(length(P),length(X));
bes = besseli(N,X);

for eP = 10.^P
    A = [A; trunc(bes)];
end

% Track number of iterations needed to achieve some accuracy

targetAccuracyAchieved = zeros(length(P),length(X));

% Try different numbers of iterations, starting with low values, and
% see how many are necessary to get desired precision

maxStart = N + 50;      % Maximum number of iterations is 50
start = N + 1;          % Which term to start on; init to single-step

epsilons = zeros(1,50); % Relative error for each number of steps

while start <= maxStart % Loop through and see when P precision is hit
    
    myi = zeros(maxStart+1,length(X));  % Clean up result array
    myi(start+1,:) = ultimateTerm;
    myi(start,:) = penultimateTerm;
    for n=start:-1:2
        for j=1:length(X)               % For each term in X...
            myi(n-1,j) = myi(n+1,j) + 2*n*myi(n,j)/X(j);
        end                             % ...recursively find In-1(x)
    end
    myi0 = myi(2,:) + 2 * myi(1,:) ./ X; % Contains I0(x) for all x in X
    
    % Normalize using 1 = I0 - 2*I2 + 2*I4 - ...
    s = myi0 * 0.5;
    for n=2:2:40
        s = s + (-1)^(n/2) * myi(n,:);
    end
    s = s .* 2;
    % Divide by s for each value in X
    for j=1:length(x)
        myi(:,j) =  nmyi ./ s;
    end
    
    % Check whether this number of iterations provided sufficient
    % decimal precision for each term. Once all terms have reached desired
    % decimal precision, skip ahead to the maxStep iteration.
    
    skipAhead = true;
    
    for i=1:length(targetAccuracyAchieved(:))
        if targetAccuracyAchieved(i) == 0
            skipAhead = false;      % Precision not reached, don't skip
            if A(i) == trunc(myi(N))
                targetAccuracyAchieved(i) = start - N - 1;
            end
        end
    end
    
    if skipAhead
        start = maxStart;       % Done honing, go for max accuracy
    else
        start = start + 1;      % Not done honing precision; do next iter
    end
    
end

% Calculate epsilon at each recursion step

for 
    
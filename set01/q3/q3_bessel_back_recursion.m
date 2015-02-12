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

% Try different numbers of iterations, starting with low values, and
% see how many are necessary to get desired precision

maxStart = N + 50;      % Maximum number of iterations is 50
start = N + 1;          % Which term to start on; init to single-step
myi = zeros(maxStart+1,length(X));

% Get actual bessel function values for besseli(N,X)

bes = zeros(maxStart + 1,length(X));
for n=1:maxStart+1
    bes(n,:) = besseli(n,X);
end

% Specify target accuracy; row indicates decimal accuracy of actual bessel
% function, column indicates x value. N is specified above.

A = zeros(length(P),length(X));

for i=1:length(P)
    for j=1:length(X)
        A(i,j) = trunc(bes(N,j),10^P(i));
    end
end

disp('Calculated actual bessel values');

% Track number of iterations needed to achieve some accuracy

targetAccuracyAchieved = zeros(length(P),length(X));

relEpsilons = zeros(maxStart+1,length(X)); % Rel error for each number of steps

while start <= maxStart % Loop through and see when P precision is hit
    
    myi = zeros(maxStart+1,length(X));  % Clean up result array
    myi(start+1,:) = ultimateTerm;
    myi(start,:) = penultimateTerm;
    for n=start:-1:2
        for j=1:length(X)               % For each term in X...
            myi(n-1,j) = myi(n+1,j) + 2*n*myi(n,j)/X(j);
        end                             % ...recursively find In-1(x)
    end
    myi0 = myi(2,:) + 2 * myi(1,:) ./ X; % Find I0(x) for all x in X
    
    % Normalize using 1 = I0 - 2*I2 + 2*I4 - ...
    s = myi0 * 0.5;         % Each x in X has it's own norm. factor
    for n=2:2:40
        s = s + (-1)^(n/2) * myi(n,:);
    end
    s = s .* 2;
    % Divide by normalization factor for each set of In(x)
    for j=1:length(X)
        myi(:,j) =  myi(:,j) / s(j);
    end
    
    % Check whether this number of iterations provided sufficient
    % decimal precision for each term. Once all terms have reached desired
    % decimal precision, skip ahead to the maxStep iteration.
    
    skipAhead = true;
    
    for p=1:length(length(P))
        for x=1:length(X)
        if targetAccuracyAchieved(p,x) == 0
            skipAhead = false;      % Precision not reached, don't skip
            if A(p,x) == trunc(myi(N,x),10^P(p))
                targetAccuracyAchieved(i) = start - N - 1;
            end
        end
        end
    end
    
    if skipAhead
        start = maxStart;       % Done honing, go for max accuracy
        disp('Done honing, go for max accuracy');
    else
        start = start + 1;      % Not done honing precision; do next iter
    end
    
end

% Calculate epsilon at each recursion step for each X

for i=1:length(X)
    relEpsilons(:,i) = abs(myi(:,i)./bes(:,i) - 1);
end

logRelEpsilons = log10(relEpsilons);

h=figure;
mesh(logRelEpsilons);
hold all;
mesh(log10(eps)*ones(size(logRelEpsilons)));
title( { 'Relative Error vs. Machine Epsilon';...
    ['While Calculating besseli(',...
    num2str(N), ',x)'] }, 'Fontsize', 18);
xlabel('x', 'Fontsize', 18);
ylabel('Steps Remaining, n', 'Fontsize', 18);
zlabel('Log10 of Relative Error', 'Fontsize', 18);
% legend({'z = log(relativeError)'; 'z = log(machineEpsilon)'}); 
print(h, '-dpdf', [ 'Error while Calculating besseli' num2str(N) '.pdf']);
function [] = cos_recursion_stefan_assignment_1()

% A simple function that recursively calculates cos(N*x)
% and plots the result vs. N. It also plots the error factor
% Nsin(N*x)/sin(x) as well as the ratio between the two.

close all
clear all

% x is the starting value for the recursion.
% Values in N are the number of steps in the recursion.
% cos(N * x) is the final value (approximately) computed.

x = 1e-6;
N = 10.^[2:7];

% Define the iteration operator Mx, which takes us from the
% vector:
%
% [cos(M * x) ; cos((M-1) * x)]
%
% to the vector:
%
% [cos((M+1) * x) ; cos(M * x)]
%
% and which is based on the cosine recursion formula:
%
% cos((M+1) * x) = 2 * cos(x) * cos(M * x) - cos((M-1) * x)

Mx = [2 * cos(x) , -1 ; 1 , 0];

% Perform the calculation for each value n in N, then store
% the result in an array of length(N) called recursionResults

recursionResults = zeros(1,length(N));

% Store the error factor n * sin(n * x) / sin(x) for each n in N 

errorFactors = zeros(1,length(N));

k=1;
for n = N
    
    % Define a vector cx to iterate over; since we are not interested
    % in intermediate steps at the moment, save resources by making
    % it a two-vector with initial value:
    %
    % [cos(x) ; cos(0)] 
    %
    % and iterate N-1 times to get final value:
    %
    % [cos(N * x) ; cos((N-1) * x)]

    cx = [cos(x) ; 1];

    % Now, calculate the above result as Mx^(N-1).
    %
    % Since I am not sure how MATLAB implements matrix powers and want
    % to avoid any funny business that would ruin my precious epsilon
    % measurements, I will implement this matrix product with a loop.

    i=1;
    while i < n                    % i and M start at 1...
        cx = Mx * cx;              % ...iterate M by +1...
        i = i + 1;                 % ...iterate i by +1...
    end                            % ...i and M end at n
    
    recursionResults(k) = cx(1);   % put cos(n * x) in result vector
    errorFactors(k) = n * sin(n * x) / sin(x);
    k = k + 1;
    
end

% Calculate ratio between recursionResults and errorFactors.

ratios = recursionResults ./ errorFactors;

h = figure('DefaultLineLineWidth', 2,'DefaultLineMarkerSize',15,...
    'DefaultTextFontSize', 18,'DefaultTextFontWeight','bold',...
    'DefaultAxesFontSize', 16);
semilogy(N,recursionResults);
hold all;
semilogy(N,errorFactors);
hold all;
semilogy(N,ratios);

% Legend string

lg = {['recursion results'];['error factor'];['ratio']};

xlabel('Number of iterations, N','FontSize',18);
ylabel('Calculated value of cos(nx)');
legend(lg,'Location','NorthWest');
hold all;

print(h, '-dpdf', [ 'cos_recursion_stefan' '.pdf' ]);
% Stefan Countryman

close all
clear all

% Run bash script to generate data

!./q2_gather_data.sh

nops = 1e4;

% Data is printed in CSV format, with each row
% indicating a new trial. Columns correspond to which
% operation is being tested.

timePerOpO0 = csvread('timing-results-o0.csv') / nops;
timePerOpO3 = csvread('timing-results-o3.csv') / nops;

% Labels for function that we tested

fxnLabels = {'',...
    'c[i] = a[i] * b[i]', 'c[i] += a[i] * b[i]', 'c[i] += a[i] * b[2*i + 20]',...
    'c[i] = a[i] / b[i]', 'c[i] += a[i] / b[i]', 'c[i] = sin(a[i])',...
    'c[i] = exp(a[i])',   'c[i] = sqrt(a[i])',   ''...
            };

% Calculate average times and plot

sizeO0 = size(timePerOpO0);
sizeO3 = size(timePerOpO3);

numTrialsO0 = sizeO0(1); % number of rows = number of trials
numTrialsO3 = sizeO3(1); % number of rows = number of trials

nFx = length(fxnLabels) - 2;
avgTimePerOpO0 = zeros(1,nFx);
avgTimePerOpO3 = zeros(1,nFx);

% Calculate avg and plot -O0 data

h = figure;
title('Compute Time for Single Operations, Optimization Level 0');
grid on; hold on; box on;
xlim([0 9]);

for i=1:nFx
    scatter(i.*ones(1,numTrialsO0),timePerOpO0(:,i),'r+');
    avgTimePerOpO0(i) = sum(timePerOpO0(:,i)) / numTrialsO0;
end

scatter(1:nFx,avgTimePerOpO0,'bx');
xlabel('Operation', 'Fontsize', 18);
ylabel('Time Per Operation (ns)', 'Fontsize', 18);
set(gca,'XTickLabel',fxnLabels);
xticklabel_rotate([],45,[],'Fontsize',12);
print(h, '-dpdf', ['Compute Time for Single Operations, GCC ' ...
                   'Optimization Level 0.pdf']);


% Calculate avg and plot -O3 data

h = figure;
title('Compute Time for Single Operations, Optimization Level 3');
grid on; hold on; box on;
xlim([0 9]);

for i=1:nFx
    scatter(i.*ones(1,numTrialsO3),timePerOpO3(:,i),'r+');
    avgTimePerOpO3(i) = sum(timePerOpO3(:,i)) / numTrialsO3;
end

scatter(1:nFx,avgTimePerOpO3,'bx');
xlabel('Operation', 'Fontsize', 18);
ylabel('Time Per Operation (ns)', 'Fontsize', 18);
set(gca,'XTickLabel',fxnLabels);
xticklabel_rotate([],45,[],'Fontsize',12);
print(h, '-dpdf', ['Compute Time for Single Operations, GCC ' ...
                   'Optimization Level 3.pdf']);

stepsInAYear = floor(365 * 24 * 3600 / r.dt / 2);

% Give each graph the same window
window = [ floor(min(min(r.y(1,:,:)))) ceil(max(max(r.y(1,:,:)))) ...
    floor(min(min(r.y(2,:,:)))) ceil(max(max(r.y(2,:,:)))) ...
    floor(min(min(r.y(3,:,:)))) ceil(max(max(r.y(3,:,:)))) ];

% Graph multiple views of the orbits to better visualize time evolution
for j=[stepsInAYear:stepsInAYear:numIterations-1 numIterations]
    h = figure('DefaultLineLineWidth',2, ...
        'DefaultTextFontSize', 18,'DefaultTextFontWeight','bold', ...
        'DefaultAxesFontSize', 16);

    formatSpec = '%7.2e';

    lg = {};    % String for the legend

    for i=1:numBodies
        positions = zeros(3,j);
        positions(:,:) = r.y(1:3,i,1:j);

        lg = vertcat(lg, ['Planet ' num2str(i) ', with mass ' ...
            num2str(mass(i)) ' solar masses']);

        plot(positions(1,:),positions(2,:));
%         plot(positions(1,:),positions(3,:));
%         plot3(positions(1,:),positions(2,:),positions(3,:));
        hold all;
    end
    axis(window(1:4));

    title({['Trajectories of ' num2str(numBodies) ' Bodies in a ' ...
        'Gravitational Field over ']; [num2str(j, formatSpec) ...
        ' Iterations Spanning ' num2str(j * r.dt / ...
        (365*24*3600)) ' years']});
    legend(lg, 'Location', 'NorthWest');
    xlabel('x');
    ylabel('y');
    % zlabel('z');

    % Print out PDF file
    print(h, '-dpdf', ['trajectories-of-' num2str(numBodies) '-bodies-in' ...
        '-a-gravitational-field-over-' num2str(j) ...
        '-iterations.pdf']);
end

h = figure('DefaultLineLineWidth',2, ...
    'DefaultTextFontSize', 18,'DefaultTextFontWeight','bold', ...
    'DefaultAxesFontSize', 16);


for energy=[ r.pet' r.ket' r.e' ]
    plot( 1:numIterations, energy );
    hold all;
end


% String for the legend
lg = {'Potential Energy'; 'Kinetic Energy'; 'Total Energy'};

title('Energy against Number of Iterations');
legend(lg, 'Location', 'NorthWest');
xlabel('Number of iterations');
ylabel('Energy');

% I've been having some trouble with the units, but one can confirm
% visually that the shapes coincide between PE & KE by looking at the
% included figures.
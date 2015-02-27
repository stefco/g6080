stepsInAYear = floor(365 * 24 * 3600 / r.dt);
for j=[stepsInAYear:stepsInAYear:numIterations-1 numIterations]
    h = figure('DefaultLineLineWidth',2, ...
        'DefaultTextFontSize', 18,'DefaultTextFontWeight','bold', ...
        'DefaultAxesFontSize', 16);

    formatSpec = '%7.2e';

    lg = {};    % String for the legend

    % I'd like to add a plot matrix with multiple
    % times and perspectives
    for i=1:numBodies
        positions = zeros(3,j);
        positions(:,:) = r.y(1:3,i,1:j);

        lg = vertcat(lg, ['Planet ' num2str(i) ', with mass ' ...
            num2str(mass(i)) ' solar masses']);

        plot3(positions(1,:),positions(2,:),positions(3,:));
        hold all;
    end

    title({['Trajectories of ' num2str(numBodies) ' Bodies in a ' ...
        'Gravitational Field over ']; [num2str(j, formatSpec) ...
        ' Iterations Spanning ' num2str(j * r.dt / ...
        (365*24*3600)) ' years']});
    legend(lg, 'Location', 'NorthWest');
    xlabel('x');
    ylabel('y');
    zlabel('z');

    % Print out PDF file
    print(h, '-dpdf', ['trajectories-of-' num2str(numBodies) '-bodies-in' ...
        '-a-gravitational-field-over-' num2str(j) ...
        '-iterations.pdf']);
end

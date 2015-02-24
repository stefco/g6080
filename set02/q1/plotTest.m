h = figure();

formatSpec = '%7.2e';

lg = {};    % String for the legend

% I'd like to add a plot matrix with multiple
% times and perspectives
for i=1:numBodies
    positions = zeros(3,numIterations);
    positions(:,:) = r.y(1:3,i,:);
    
    lg = vertcat(lg, ['Planet ' num2str(i) ', with mass ' ...
        num2str(mass(i)) ' solar masses']);
    
    plot3(positions(1,:),positions(2,:),positions(3,:));
    hold all;
end

title(['Trajectories of ' num2str(numBodies) ' Bodies in a ' ...
    'Gravitational Field over ' num2str(numIterations) ...
    ' Iterations Spanning ' num2str(r.dt) ' seconds']);
legend(lg, 'Location', 'NorthWest');
xlabel('x');
ylabel('y');
zlabel('z');
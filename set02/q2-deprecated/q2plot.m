h = figure('DefaultLineLineWidth',2, ...
        'DefaultTextFontSize', 18,'DefaultTextFontWeight','bold', ...
        'DefaultAxesFontSize', 16);

    formatSpec = '%7.2e';

    lg = {};    % String for the legend
    
    mass  = 39.948 / 6.0221e23 * 1e-3;
    
    % Grab particle velocities from every 10 steps and plot them
    [ np, numberOfKEMeasurements ] = size( r.ke );
    sqVelocities = [];
    for i=50:10:numberOfKEMeasurements  % Give a bit of time to thermalize
        sqVelocities = [ sqVelocities; r.ke( :, i ) ];
    end
    
    % Take square roots
    velocities = sqVelocities .^(1/2);
        
    histogram(velocities, 20)
    hold all;
    
    % Pick velocity range that overlaps with simulation results
    v = 0:0.001:1 * max( velocities );
    
    % Calculate Maxwell distribution; can't get it plotting properly on
    % top of the histogram
    f = ( 3 * mass / 4 * pi * r.ket(end) ) ^ ( 3 / 2 ) ...
            * 4 * pi .* (v.^2) .* exp( -3 * mass * (v.^2) ...
            / ( 4 * r.ket(end) ) );

    plot( v, f )

    title({['Velocity distribution of Argon particles vs. ' ...
        'Maxwell Distribution with Same Parameters']});
%     legend(lg, 'Location', 'NorthWest');
    xlabel('t');
    ylabel('density');

    % Print out PDF file
    print(h, '-dpdf', ['velocity-distribution-of-argon-particles-vs-' ...
        'maxwell-distribution.pdf']);

h = figure('DefaultLineLineWidth',2, ...
    'DefaultTextFontSize', 18,'DefaultTextFontWeight','bold', ...
    'DefaultAxesFontSize', 16);


for energy=[ r.pet(1:100)' r.ket(1:100)' r.e(1:100)' ]
    plot( 1:100, energy );
    hold all;
end

% String for the legend
lg = {'Potential Energy'; 'Kinetic Energy'; 'Total Energy'};

title('Energy against Number of Iterations');
legend(lg, 'Location', 'NorthWest');
xlabel('Number of iterations');
ylabel('Energy');
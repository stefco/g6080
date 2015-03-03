% Part 3
numBodies       = 3;

% (Again), my energy tracking is off; I had trouble figuring out
% the units. The graphs have the same forms as one another, so I know
% it's just a scaling issue. That said, there is very clearly a completely
% unphysical, extremely sudden spike in kinetic energy as a result of the
% too-short timescale.

% However, we can address the problem conceptually using analytical
% reasoning. Since the acceleration scales like 1/r^2, the maximum is much
% larger for this trial (due to the close passes). This means that a Taylor
% series in time, with a now-huge acceleration factor for the t^2 power,
% would no longer be linear for what might have previously been a
% sufficiently small time scale. Without resorting to higher-order
% approximations, the only choice for maintaining accuracy is to use a
% smaller time step.

% Compare the results of Part 3a.
numIterations   = 2.5e4;
X               = 2e-4;
mass            = [ 1 1e-2 1e-2 ];
initial         = zeros( 6, numBodies );
initial( 1, 2 ) = 1;
initial( 5, 2 ) = 1;
initial( 1, 3 ) = -1;
initial( 5, 3 ) = 1.01;

r = planetaryMotion( numBodies, numIterations, X, mass, initial );

plotTest;
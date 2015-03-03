% Part 3 (continued)
numBodies       = 3;

% Smaller timestep leads to much smoother behaviour; decreasing it by an
% order of magnitude removes the most egregiously anomalous energy
% conservation violations and produces a smooth, albeit tight, bend where
% the planets "slingshot" off of one another.
numIterations   = 2.5e5;
X               = 2e-5;
mass            = [ 1 1e-2 1e-2 ];
initial         = zeros( 6, numBodies );
initial( 1, 2 ) = 1;
initial( 5, 2 ) = 1;
initial( 1, 3 ) = -1;
initial( 5, 3 ) = 1.01;

r = planetaryMotion( numBodies, numIterations, X, mass, initial );

plotTest;
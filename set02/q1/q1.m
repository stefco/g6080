% Part 1
numBodies       = 3;
numIterations   = 9e4;
X               = 2e-4;
mass            = [ 1 1e-2 1e-2 ];
initial         = zeros( 6, numBodies );
initial( 1, 2 ) = 1;
initial( 5, 2 ) = 1;
initial( 1, 3 ) = -1;
initial( 5, 3 ) = 1.04;

r = planetaryMotion( numBodies, numIterations, X, mass, initial );

plotTest;
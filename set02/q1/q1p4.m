% Part 4
numBodies       = 4;

% After a bunch of tweaking, this was one of my better "slingshots";
% suffice it to say, it's a hard problem to eyeball. Better to write a 
% program to iteratively seek out stronger slingshot boosts.
numIterations   = 0.2e5;
X               = 2e-5;
mass            = [ 1 1e-2 0.5e-2 1e-5 ];
initial         = zeros( 6, numBodies );
initial( 1, 2 ) = 1;
initial( 5, 2 ) = 1;
initial( 1, 3 ) = 2/3;
initial( 5, 3 ) = (2/3)^(-0.5);
initial( 1, 4 ) = 0.99;
initial( 4, 4 ) = -2.07;
initial( 5, 4 ) = 1.19;

r = planetaryMotion( numBodies, numIterations, X, mass, initial );

plotTest;
% Part 2
numBodies       = 3;

% Within 18,000 iterations, we see a slingshot effect; the outer
% planet gets sufficiently close to the inner planet to "fling" it out into
% a higher orbit.

% Due to momentum conservation (and the fact that they provide the dominant
% forces acting on one another during collision), the planets end up flying
% off in near-opposite directions; this puts them on track to collide
% again, since they're both in nearly the same elliptical orbit. This
% result relies on their equal masses to turn momentum conservation into
% equal/opposite velocities.

% Unfortunately, my energy tracking is off; I had trouble figuring out
% the units. The graphs have the same forms as one another, so I know
% it's just a dumb scaling issue.
numIterations   = 1.8e4;
X               = 2e-4;
mass            = [ 1 1e-2 1e-2 ];
initial         = zeros( 6, numBodies );
initial( 1, 2 ) = 1;
initial( 5, 2 ) = 1;
initial( 1, 3 ) = -1;
initial( 5, 3 ) = 1.01;

r = planetaryMotion( numBodies, numIterations, X, mass, initial );

plotTest;
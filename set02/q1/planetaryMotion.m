% 2/22/2015 Stefan Countryman

function [ r ] = planetaryMotion( numBodies, numIterations, X, mass, initial )
% PLANETARYMOTION  1st-order planetary motion predictor-corrector
%   r = PLANETARYMOTION(p) takes a struct p where, e.g.,
%
%     NB             (number of planets/bodies)
%     N              (number of simulation steps)
%     X              (step size)
%     mass           (masses in solar mass units);f
%     init(d,i)      (init values, dth dimension, ith body,
%                      rEarth units)
%
%   And returns a struct r with those same properties as well as
%
%     r.y(d,i,n)     (values, dth dimension, ith body, nth step)
%     r.f(d,i,n)     (forces, dth dimension, ith body, nth step)
%     r.pe(i,n)      (potential energy,      ith body, nth step)
%     r.ke(i,n)      (kinetic energy,        ith body, nth step)
%     r.pet(n)       (total potential energy,          nth step)
%     r.ket(n)       (total kinetic energy,            nth step)
%     r.e(n)         (total energy,                    nth step)

% Close all open graphs
close all

% Set constants
rEarth = 1.495e11;     % Earth mean orbital radius, m
mSolar = 1.9891e30;    % Solar mass, kg
mEarth = 5.972e24;     % Earth mass, kg
gNewt  = 6.673e-11;    % Gravitational Constant, N m^2 kg^-2

% Determine natural units
tUnits = sqrt( rEarth^3 / ( gNewt * mSolar ) );
gTilde = gNewt * mSolar * tUnits^2 / rEarth^4;

% Check arguments
if (numBodies < 2)
    error('planetaryMotion:numIterations', ['numBodies must be an ' ...
                        'int >= 2.']);
elseif (numIterations < 1)
    error('planetaryMotion:numIterations', ['numIterations must be a ' ...
                        'positive int.']);
elseif (X < 0)
    error('planetaryMotion:stepSize', 'X must be positive.');
elseif (size(mass) ~= [1 numBodies])
    error('planetaryMotion:argumentMismatch', ['mass must be ' ...
                        'row vector with length numBodies.']);
elseif (size(initial) ~= [6 numBodies])
    error('planetaryMotion:argumentMismatch', ['init must ' ...
                        'have 6 rows and numBodies columns.']);
end

% Make sure masses are positive
arePositive = (mass > 0);
for isPositive=arePositive
    if ~isPositive
        error('planetaryMotion:negativeMasses', ['mass ' ...
                            'values must be positive']);
    end
end

% Initialize results
r.NB  = numBodies;            % Number of bodies, >=2
r.N   = numIterations;        % Number of simulation steps
r.dx  = X;                    % Step size
r.dt  = X * tUnits;           % Step size, in seconds
r.m   = mass;                 % Masses of bodies

r.y   = zeros(6,numBodies, numIterations);
r.f   = zeros(6,numBodies, numIterations);
r.pe  = zeros(numBodies, numIterations);
r.ke  = zeros(numBodies, numIterations);
r.pet = zeros(1, numIterations);
r.ket = zeros(1, numIterations);
r.e   = zeros(1, numIterations);
r.y(:,:,1) = initial;

% Forces and potential/kinetic/total energies for 1st step
KE(1);
PEandF(1);
r.e(1) = r.ket(1) + r.pet(1);

% Loop through all iterations
for n = 1:numIterations-1
  % Predict position:    Yp[n+1] = Y[n] + X*F[n]
  r.y(:,:,n+1) = r.y(:,:,n) + X * r.f(:,:,n);
  % Predict derivative:  Fp[n+1] = F(Yp(n+1))
  PEandF(n+1);
  % Correct position:    Y[n+1]  = Y[n] + X*Fp[n+1]
  r.y(:,:,n+1) = r.y(:,:,n) + X * r.f(:,:,n+1);
  % Correct derivative:  F[n+1]  = F(Y(n+1))
  PEandF(n+1);
  % Kinnetic energy
  KE(n+1);
  % Total energy
  r.e(n+1) = r.ket(n+1) + r.pet(n+1);
end

% Save results file
save('planetary-motion-results', 'r');

% Function for calculating Kinetic Energy
function [] = KE(nn)

  % nnth total kinetic energy starts at 0 
  r.ket(nn) = 0.0;

  % Calculate ke for each body
  for ii=1:numBodies
    r.ke(ii,nn) = 0.5 * r.m(ii) ...
      * dot(r.y(4:6,ii,nn),r.y(4:6,ii,nn));
    r.ket(nn) = r.ket(nn) + r.ke(ii,nn);
  end
end

% Function for calculating Potential Energy and
% derivatives of each variable
function [] = PEandF(nn)

  % nnth total potential energy starts at zero
  r.pet(nn) = 0.0;

  % Calculate pe and f for each body
  for ii=1:numBodies
    % Initialize to zero
    r.pe(ii,nn) = 0;
    r.f(4:6,ii,nn) = zeros(3,1);

    % Get dr/dx components for F[nn]
    r.f(1:3,ii,nn) = r.y(4:6,ii,nn);
    
    % Get contributions from other bodies
    for jj=1:numBodies
      if jj ~= ii
        % Calculate 1/r^2
        deltaR = r.y(1:3,ii,nn) - r.y(1:3,jj,nn);
        normR = norm( deltaR );
        % Add jjth potential energy contribution
        r.pe(ii,nn) = r.pe(ii,nn) - gTilde ...
            * mass(ii) * mass(jj) / normR;
        % Add jjth dv/dx contribution
        r.f(4:6,ii,nn) = r.f(4:6,ii,nn) - mass(jj) ...
            / normR^3 .* deltaR;
      end
    end

    % Add energy contribution to total PE
    r.pet(nn) = r.pet(nn) + r.pe(ii,nn);
  end
end
end

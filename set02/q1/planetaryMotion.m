% 2/22/2015 Stefan Countryman

function [ r ] = planetaryMotion( varargin )
% PLANETARYMOTION  1st-order planetary motion predictor-corrector
%   r = PLANETARYMOTION(p) takes a struct p where, e.g.,
%
%     p.NB           (number of planets/bodies)
%     p.N            (number of simulation steps)
%     p.dx           (step size)
%     p.m            (masses in solar mass units)
%     p.in(d,i)      (init values, dth dimension, ith body,
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
gTilde = G * mSolar^2 * rEarth^2;
tUnits = sqrt( rEarth^3 / ( gNewt * mSolar ) );

if nargin = 0
    % Initialize default values, from RDM's assignment
    p.NB = 3;
    p.N  = 3e4;
    p.dx = 2e-4;
    p.m  = [1, 1e-2, 1e-2];
    % Initial values for bodies (note that sun rests at center):
    p.in = zeros(6, p.NB);
    p.in(1,2) = 1;
    p.in(5,2) = 1;
    p.in(1,3) = -1;
    p.in(5,3) = 1.04;
end

% Check arguments
if ~isinteger(p.NB) | p.NB < 2
    error('planetaryMotion:numberBodies', ['p.NB must be an ' ...
                        'int >= 2.']);
elseif ~isinteger(p.N) | p.NB < 1
    error('planetaryMotion:numberSteps', ['p.NB must be a ' ...
                        'positive in.']);
elseif p.dx < 0
    error('planetaryMotion:stepSize', 'p.dx must be positive.');
elseif size(p.m) ~= [1 p.NB]
    error('planetaryMotion:argumentMismatch', ['p.m must be ' ...
                        'row vector with length p.NB.']);
elseif size(p.in) ~= [6 p.NB]
    error('planetaryMotion:argumentMismatch', ['p.in must ' ...
                        'have 6 rows and p.NB columns.']);
end

% Make sure masses are positive
arePositive = (p.m > 0);
for isPositive=arePositive
    if ~isPositive
        error('planetaryMotion:negativeMasses', ['p.m ' ...
                            'values must be positive']);
    end
end

% Use local variables for simplicity
numBodies = p.NB;
numIterations = p.N;
X = p.dx;
mass = p.m;
r.y(:,:,1) = p.in;

% Initialize results
r.NB  = p.NB;                 % Number of bodies, >=2
r.N   = p.N;                  % Number of simulation steps
r.dx  = p.dx;                 % Step size
r.m   = p.m;                  % Masses of bodies

r.y   = zeros(6,numBodies, numIterations);
r.f   = zeros(6,numBodies, numIterations);
r.pe  = zeros(numBodies, numIterations);
r.ke  = zeros(numBodies, numIterations);
r.pet = zeros(1, numIterations);
r.ket = zeros(1, numIterations);
r.y   = zeros(1, numIterations);

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
  r.y[:,:,n+1] = r.y[:,:,n] + X * r.f[:,:,n+1];
  % Correct derivative:  F[n+1]  = F(Y(n+1))
  PEandF(n+1);
  % Kinnetic energy
  KE(n+1);
  % Total energy
  r.e(n+1) = r.ket(n+1) + r.pet(n+1);
end

% Function for calculating Kinetic Energy
function [] = KE(nn)

  % nnth total kinetic energy starts at 0 
  r.ket(nn) = 0.0;

  % Calculate ket for each body
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

  % Declare variables
  deltaR = zeros(3,1);
  normR = 0.0;
  
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

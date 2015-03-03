% 2/25/2015

function [ r ] = molecularDynamics( numBodies, particleDensity, ...
    temperature, numIterations, varargin )

% Close all open graphs
close all

nargin

% Set constants
sigma = 3.405e-10;   % Length units, from Verlet, in m
epsK  = 119.8;       % Energy units, from Verlet, in Kelvin
kB    = 1.3806e-23;  % Boltzmann constant, in J / K
eps   = 1.5 * kB * epsK; % Energy units, from Verlet, in Joules
mass  = 39.948 / 6.0221e23 * 1e-3; % Argon mass, in kg
tau   = sqrt( mass * sigma^2 / ( 48 * eps ) ); % Time units, in s
indexInterval = 10;  % How often to index particles for force calc
saveInterval  = 10;  % How often to save positions, velocities, etc.
numGridDivs   = 10;  % Lengthwise subdivisions of box for indexing
includeDivs   = 2;   % How many adjacent subdivisions to use for force
sqInterDist   = 4;   % Square of maximum worthwhile interaction distance

% Make it possible to resume from where we left off by loading a file
if nargin == 5
    resultsFileName = varargin{1};
else

    % Set start date
    startDate = datetime('now');
    startDate.Format = 'yyyy-MM-dd-HH-mm-ss';
    resultsFileName = [ char(startDate) '-trial' ];


    % Initialize results variable
    r.numBodies  = numBodies;             % Number of particles
    r.particleDensity = particleDensity;  % Density of particles
    r.temperature = temperature;          % Desired temperature, in K
    r.dx  = 0.032;                        % Time step used by Verlet, in tau
    r.dt  = r.dx * tau;                   % Time step used by Verlet, in s
    r.y   = zeros( 3, numBodies, 1 );     % Positions
    r.v   = zeros( 3, numBodies, 1 );     % Velocities
    r.f   = zeros( 3, numBodies, 1 );     % Accelerations
    r.pe  = zeros( numBodies, 1 );        % Potential Energies
    r.ke  = zeros( numBodies, 1 );        % Kinetic Energies
    r.pet = zeros( 1 );                   % Total Potential Energy
    r.ket = zeros( 1 );                   % Total Kinetic Energy
    r.e   = zeros( 1 );                   % Total Energy
    r.T   = zeros( 1 );                   % Temperature

    r.volume = 0;
    r.sideLength = 0;
    r.saveCount = 1;
    r.indexCount = 1;

    r.subDivs = zeros( 3, numBodies );    % Sort particles into boxes
    r.index = cell( numGridDivs, numGridDivs, numGridDivs );
    r.interactionList = [];

    % Place particles in a uniform grid, with slight random fluctuations
    mdPlace();
end
    
disp(['On Step 1 of ' num2str(numIterations)]);

% Index particles for performance boost
index(1);

% Calculate energies for 1st step
KE( 1 );
PEandF(1);
r.e(1) = r.ket(1) + r.pet(1);

% Make the temperature equal the desired temperature
% mdAdjustTemp();

% Iterate for the desired number of iterations
for nn = 1:numIterations-1
    disp(['On Step ' num2str(nn + 1) ' of ' num2str(numIterations)]);
    % Tack on matrix entries for positions and velocities
    r.y = cat( 3, r.y, zeros( 3, numBodies ) );
    r.v = cat( 3, r.v, zeros( 3, numBodies ) );
    % Calculate next position
    r.y( :, :, nn + 1 ) = r.y( :, :, nn ) + r.v( :, :, nn ) ...
        * r.dx + 0.5 * r.dx^2 * r.f( :, :, nn );
    % If position is past boundaries, bring back to other side
    r.y( :, :, nn + 1 ) = mod( r.y( :, :, nn + 1 ), r.sideLength );
    % Calculate next potential energy and acceleration
    PEandF( nn + 1 );
    % Calculate intermediate velocity
    r.v( :, :, nn + 1 ) = r.v( :, :, nn ) + 0.5 * r.dx ... 
        * r.f( :, :, nn );
    % Calculate next velocity
    r.v( :, :, nn + 1 ) = r.v( :, :, nn + 1 ) + 0.5 * r.dx ...
        * r.f( :, :, nn );
    % Calculate next kinetic energy
    KE( nn + 1 );
    % Total energy
    r.e = [ r.e ( r.ket( nn + 1 ) + r.pet( nn + 1 ) ) ];
    % See if it's time to save
    if r.saveCount < saveInterval
        r.saveCount = r.saveCount + 1;
    else
        % This struct contains the entire state of the simulation at this
        % timestep; we'll use its contents later to calculate pressure.
        save( [ resultsFileName '-r' ], 'r' );
        r.saveCount = 1;
    end
    % See if it's time to index
    if r.indexCount < indexInterval
        r.indexCount = r.indexCount + 1;
    else
        index( nn + 1 );
        r.indexCount = 1;
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%  AUXILIARY FUNCTIONS
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    % FUNCTION FOR INITIALLY PLACING PARTICLES
	function [ ] = mdPlace()
	
        % Find volume and, subsequently, side lengths for box
        r.volume = r.numBodies / r.particleDensity;
        r.sideLength = r.volume ^ (1/3);

        % Fit particles to evenly spaced grid. There might be gaps.
        gridNotches = ceil( r.numBodies ^ (1/3) );
        gridUnit = r.sideLength / gridNotches;

        for ii = 1:r.numBodies
            num = ii - 1; % Start at zero; none of this MATLAB 1 nonsense
            xGrid = mod( num, gridNotches );
            num = num - xGrid;
            yGrid = mod( num, gridNotches^2 ) / gridNotches;
            num = num - yGrid * gridNotches;
            zGrid = num / gridNotches^2;
            r.y( :, ii, 1 ) = gridUnit .* [ xGrid; yGrid; zGrid ];
        end

        % Place particles in middle of each cell, with some random variation
        r.y( :, :, 1 ) = r.y( :, :, 1 ) + gridUnit .* ...
                ( 0.1 * rand( 3, r.numBodies ) + 0.45 );

        disp( 'done placing particles' );
        fname = [ char( startDate ) '-place' ];
        dlmwrite( fname, r.y);
        disp( [ 'wrote particle placements to ' fname ] );

    end
	
    % FUNCTION FOR ORGANIZING PARTICLES TO REDUCE FORCE CALCULATION TIME
    function [ ] = index( nnn )
        
        disp( [ 'Generating index for step ' num2str( nnn ) ] );
        
        % Reset interaction list and inhabitants of each subdivision
        r.index = cell( numGridDivs, numGridDivs, numGridDivs );
        r.interactionList = [];
        
        % Assign each particle to a subdivision of the box
        r.subDivs = ceil( r.y( :, :, nnn ) * numGridDivs / r.sideLength );
        
        % Place each particle in a subdivision
        for bb=1:numBodies
            
            ys = r.subDivs( :, bb );
            r.index{ys(1),ys(2),ys(3)} ...
                = [ r.index{ys(1),ys(2),ys(3)} bb ];
        end
        
        % Create an "interaction list" which pairs particles that
        % are close enough to interact strongly, using nearby candidates
        
        YS = zeros( 3, 2 * includeDivs + 1 ); % Coords of close subdivs

        for bb=1:numBodies
            % The location of body bb
            ys = r.subDivs( :, bb );
            % All the nearby boxes
            YS(1,:) = ys(1)-includeDivs:ys(1)+includeDivs;
            YS(2,:) = ys(2)-includeDivs:ys(2)+includeDivs;
            YS(3,:) = ys(3)-includeDivs:ys(3)+includeDivs;
            YS = mod( YS, numGridDivs ) + 1;
            % Look in each box for other bodies
            for y1=YS(1,:); for y2=YS(2,:); for y3=YS(3,:)
                for ob=r.index{y1,y2,y3}
                    % No self-interaction
                    if ob == bb
                        continue
                    end
                    sqd = rSquaredAndDeltaR(...
                        r.y( :, bb, nnn )', r.y( :, ob, nnn )' );
                    if sqd < sqInterDist
                        % These particles are near, and will interact
                        r.interactionList = ...
                            [r.interactionList [bb;ob]];
                    end
                end
            end; end; end
        end
        
        disp([ 'Interaction list has ' num2str(length ...
            (r.interactionList) ) ' pairs (with duplication).'] );
        
    end

    % FUNCTION FOR GETTING SQUARE DISTANCE AND SHORTEST VECTOR BETWEEN
    % TWO POINTS
    function [sqd, del] = rSquaredAndDeltaR(point2,point1)
        
        % This lets us take col vectors as input without trouble
        p2 = zeros( 1, 3 );
        p1 = zeros( 1, 3 );
        p2(:) = point2(:);
        p1(:) = point1(:);
        
        % Calculate the minimum squared distance, accounting for
        % the periodicity of the containing box. Note that this will
        % also return the index of the min for each column.
        deltaCandidates = [ (p2-p1); (p2-(p1-r.sideLength)); ...
            (p2-(p1+r.sideLength)) ];
        [sqds, ind] = min( deltaCandidates.^2 );
        sqd = sum(sqds);
        % Each column has the "long" and "short" distance in that dimension
        del = zeros(1,3);
        for iii=1:3
            % For each column, pick the shorter-distance row
            del(iii) = deltaCandidates(ind(iii), iii);
        end
    end

    % FUNCTION FOR CALCULATING POTENTIAL ENERGY AND ACCELERATION
    function [] = PEandF( nnn )

      % nnnth potential energies and forces start at zero
      r.f   = cat( 3, r.f, zeros( 3, numBodies ) );
      r.pet = [ r.pet 0.0 ];
      r.pe  = [ r.pe zeros( numBodies, 1 ) ];

      % Calculate rSquared, pe, and f for each interacting body; recall
      % that ii(1) is the current body, and ii(2) is the other body.
      for ii=r.interactionList
        % Get rSquared and deltaR
        [rSquared, deltaR] = rSquaredAndDeltaR( r.y( :, ii(1), nnn ), ...
            r.y( :, ii(2), nnn ) );
        
        % Compute acceleration from ii(2) particle
        extraAccel = deltaR * ( 1/rSquared^7 - 0.5 * 1/rSquared^4 );
        
        % Add acceleration
        r.f( :, ii(1), nnn ) = r.f( :, ii(1), nnn ) + extraAccel(:);

        % Compute potential energy from ii(2) particle
        extraPE = 4 * eps * ( 1/rSquared^6 - 1/rSquared^3 );
        
        % Add potential energy
        r.pe( ii(1), nnn ) = r.pe( ii(1), nnn ) + extraPE(:);

        % Add energy contribution to total PE
        r.pet(nnn) = r.pet(nnn) + r.pe(ii(1),nnn);
        
      end
      
    end
	
    % FUNCTION FOR CALCULATING KINETIC ENERGY
    function [] = KE( nnn )

      % nnnth kinetic energies start at 0 
      r.ket = [ r.ket 0.0 ];
      r.ke  = [ r.ke zeros( numBodies, 1 ) ];
      r.T   = [ r.T 0.0 ];

      % Calculate ke for each body
      for ii=1:numBodies
        % Converting from our natural units to SI, we get a factor of
        % m * sigma^2 / tau^2, which equals 48 * eps by earlier requirement
        r.ke(ii,nnn) = 24 * eps ...
          * dot(r.v(:,ii,nnn),r.v(:,ii,nnn));
        r.ket(nnn) = r.ket(nnn) + r.ke(ii,nnn);
      end
      r.T(nnn) = 2 * ket(nnn) / ( 3 * r.numBodies * kB );
    end

% 	function [ ] = mdAdjustTemp()
% 	
% 	% Iterate 25 steps forward
% 	mdIterate( 25 );
%     
%     % 
% 	
% 	end
	
end
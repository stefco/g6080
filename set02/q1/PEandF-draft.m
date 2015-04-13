% 2/23/15 Stefan Countryman

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

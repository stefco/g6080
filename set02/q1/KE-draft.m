% 2/23/15 Stefan Countryman

function [] = KE(nn)

  % nnth total kinetic energy starts at 0 
  r.ket(nn) = 0.0;

  % Calculate ke for each body
  for ii=1:numBodies
    r.ke(ii,nn) = 0.5 * mass(ii) ...
      * dot(r.y(4:6,ii,nn),r.y(4:6,ii,nn));
    r.ket(nn) = r.ket(nn) + r.ke(ii,nn);
  end
end

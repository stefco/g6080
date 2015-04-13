%  2/13/15  RDM
%
%  First order predictor/corrector for planetary motion

function [ p ] = planetary_motion()

close all

p.NP = 3;       %  Number of particles/planets >= 2
p.N = 30000;    %  Number of steps
p.dx = 2e-4;            %  step size

p.m = [1 ,0.01, 0.01];    %  Mass of particles/planets

if ( size(p.m, 2) ~= p.NP )
    error('Particle number %d and number of masses %d differ\n', ...
        NP, size(p.m, 2));
end

%  Local variables for clarity
NP = p.NP;
N = p.N;

p.y = zeros(6,NP,N);    %  position 1:3 and velocity 4:6
p.f = zeros(6,NP,N);    %  forces
p.pe = zeros(NP,N);     %  potential energy for particle i
p.ke = zeros(NP,N);     %  kinetic energy for particle i
p.pet = zeros(1,N);     %  total potential energy
p.ket = zeros(1,N);     %  total kinetic energy
p.e = zeros(1,N);       %  totel energy

%  Initial conditions
p.y(1,2,1) = 1;
p.y(5,2,1) = 1;
p.y(1,3,1) = -1;
p.y(5,3,1) = 1.04;

KE(1);
PEandF(1);
p.e(1) = p.ket(1) + p.pet(1);

for n = 1:N-1
    p.y(:,:,n+1) = p.y(:,:,n) + p.dx * p.f(:,:,n);
    PEandF(n+1);
    p.y(:,:,n+1) = p.y(:,:,n) + p.dx * p.f(:,:,n+1);
    PEandF(n+1);
    KE(n+1);
    p.e(n+1) = p.ket(n+1) + p.pet(n+1);
end
   
figure
plot(squeeze(p.y(1,1,:)), squeeze(p.y(2,1,:)));
legend('Particle 1');
hold all

figure
for m=2:NP
plot(squeeze(p.y(1,m,:)), squeeze(p.y(2,m,:)));
hold all
end

figure
plot(p.pet,'b');
hold all
plot(p.ket,'g');
hold all
plot(p.pet + p.ket,'r');
legend(['PE', 'KE', 'Total E' ]);


%  Find potential energy and force for particle i

    function [] = PEandF(nn)
  
    end

    function [] = KE(nn)
        
        p.ket(nn) = 0.0;
        
        for ii=1:NP
            
            p.ke(ii,nn) = 0.5 * p.m(ii) * ...
                ( p.y(4,ii,nn) * p.y(4,ii,nn) ...
                +  p.y(5,ii,nn) * p.y(5,ii,nn) ...
                + p.y(6,ii,nn) * p.y(6,ii,nn));
            
            p.ket(nn) = p.ket(nn) + p.ke(ii,nn);
        end
        
    end

end


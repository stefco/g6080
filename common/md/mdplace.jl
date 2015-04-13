# Fit particles to evenly spaced grid. There might be unused spots.
function mdplace!(r::MolecularDynamicsTrial)
    # Make sure this isn't a trial in progress
    if r.currentStep != 0
        error("Cannot place particles: currentStep = ", r.currentStep);
    end

    # We need at least enough grid notches for the number of bodies; round up
    gridNotches = ceil( r.numBodies ^ (1/3) );

    # We need to know how long each grid unit is
    gridUnit = r.sideLength / gridNotches;

    # Iterate through the bodies and lay them down in neat little rows
    for i in 1:r.numBodies
        num = i - 1; # Because Julia and MATLAB start at 1, not 0

        # Provide jth coordinate value and shift radix representation right
        for j in 1:3 
            grid = ( num % gridNotches );       # find radix^j place
            num = ( num - grid ) / gridNotches; # shift right
            grid += 0.1rand() + 0.45;           # middle of cell, some variation
            r.y[j,i,1] = grid * gridUnit;       # scale by box size
        end
    end

    # Particles are now placed, so say we are on step 1
    r.currentStep = 1;
end
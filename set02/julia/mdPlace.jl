# Fit particles to evenly spaced grid. There might be gaps.
function mdplace(r::MolecularDynamicsTrial)
    # Make sure this isn't a trial in progress
    if r.currentStep != 0
        error("Cannot place particles: currentStep = ", r.currentStep);
    end

    # We need at least enough grid notches for the number of bodies; round up
    gridNotches = ceil( r.numBodies ^ (1/3) );

    # We need to know how long each grid unit is
    gridUnit = r.sideLength / gridNotches;

    # Iterate through the bodies and lay them down in neat little rows
    for i = 1:r.numBodies
        # Because Julia and MATLAB start at 1, not 0
        num = i - 1; 

        # Provide x position
        grid = ( num % gridNotches );
        num = ( num - grid ) / gridNotches;
        r.y[1,i,1] = grid * gridUnit;

        # Provide y position
        grid = ( num % gridNotches );
        num = ( num - grid ) / gridNotches;
        r.y[2,i,1] = grid * gridUnit;

        # Provide z position
        grid = num;
        r.y[3,i,1] = grid * gridUnit;
    end

    # Place particles in middle of each cell, with some random variation
    r.y[:,:,1] += ( 0.1 * rand( 3, r.numBodies ) + 0.45 ) * gridUnit;

    println("Particles placed.");
end

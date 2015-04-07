type MolecularDynamicsTrial
    # Initial Parameters
    numBodies::Int64                # Number of particles
    particleDensity::Float64        # Density of particles
    desiredTemperature::Float64     # Desired temperature
    dx::Float64                     # Time step; defaults to Verlet's value
    volume::Float64                 # Volume of box
    sideLength::Float64             # Side length of box
    numIterations::Int64            # How long to run the simulation

    # Dynamic Quantities
    y::Array{Float64, 3}            # Positions
    v::Array{Float64, 3}            # Velocities
    f::Array{Float64, 3}            # Accelerations
    pe::Array{Float64, 2}           # Potential Energies
    ke::Array{Float64, 2}           # Kinetic Energies
    pet::Array{Float64, 1}          # Total Potential Energy
    ket::Array{Float64, 1}          # Total Kinetic Energy
    e::Array{Float64, 1}            # Total Energy
    vir::Array{Float64, 1}          # Virial, ∑∑r(∂V/∂r)
    currentStep::Int64              # Where this simulation left off
end

# Outer constructor just for input parameters
MolecularDynamicsTrial
        (numBodies, particleDensity, desiredTemperature, dx, numIterations) 
        = MolecularDynamicsTrial(
            numBodies,
            particleDensity,
            desiredTemperature,
            dx,
            (numBodies / particleDensity),          # Volume
            (numBodies / particleDensity)^(1/3),    # Side length
            numIterations,
            zeros( 3, numBodies, numIterations ),   # Positions
            zeros( 3, numBodies, numIterations ),   # Velocities
            zeros( 3, numBodies, numIterations ),   # Accelerations
            zeros( numBodies, numIterations ),      # Potential Energies
            zeros( numBodies, numIterations ),      # Kinetic Energies
            zeros( numIterations ),                 # Total Potential Energy
            zeros( numIterations ),                 # Total Kinetic Energy
            zeros( numIterations ),                 # Total Energy
            zeros( numIterations ),                 # Virial
            0                                       # Completed steps (start at 0)
        );


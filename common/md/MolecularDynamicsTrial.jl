# Struct for holding all parameters associated with an md simulation
type MolecularDynamicsTrial
    # Parameters
    numBodies::Int64                # Number of particles
    ρ::Float64                      # Density of particles
    Td::Float64                     # Desired temperature
    dx::Float64                     # Time step; defaults to Verlet's value
    V::Float64                      # Volume of box
    sideLength::Float64             # Side length of box
    steps::Int64                    # How long to run the simulation
    σ::Float64                      # Length scale, from Verlet, in m
    ϵ::Float64                      # Energy scale, in Kelvin
    m::Float64                      # Argon atom mass, in kg
    τ::Float64                      # Time scale, in seconds

    # Dynamic Quantities
    y::Array{Float64, 3}            # Positions
    v::Array{Float64, 3}            # Velocities
    f::Array{Float64, 3}            # Accelerations
    pe::Array{Float64, 2}           # Potential Energies
    ke::Array{Float64, 2}           # Kinetic Energies
    pet::Array{Float64, 1}          # Total Potential Energy
    ket::Array{Float64, 1}          # Total Kinetic Energy
    T::Array{Float64, 1}            # Temperature at each step
    e::Array{Float64, 1}            # Total Energy
    vir::Array{Float64, 1}          # Virial, ∑∑r(∂V/∂r)
    currentStep::Int64              # Where this simulation left off
end

# Outer constructor just for input parameters
MolecularDynamicsTrial( numBodies, ρ, Td, dx, steps, σ, ϵ, m ) =
    MolecularDynamicsTrial(
        numBodies,
        ρ,
        Td,
        dx,
        ( numBodies / ρ ),                      # Volume
        ( numBodies / ρ )^( 1 / 3 ),            # Side length
        steps,
        σ,                                      # Length scale, in m
        ϵ,                                      # Energy scale, in Kelvin
        m,                                      # Argon atom mass, in kg
        sqrt( m * σ^2 / ( 48ϵ * 1.3806e-23 ) ), # Time scale, in seconds
        zeros( 3, numBodies, steps ),           # Positions
        zeros( 3, numBodies, steps ),           # Velocities
        zeros( 3, numBodies, steps ),           # Accelerations
        zeros( numBodies, steps ),              # Potential Energies
        zeros( numBodies, steps ),              # Kinetic Energies
        zeros( steps ),                         # Total Potential Energy
        zeros( steps ),                         # Total Kinetic Energy
        zeros( steps ),                         # Temperature at each step
        zeros( steps ),                         # Total Energy
        zeros( steps ),                         # Virial
        0 );                                    # Completed steps (start at 0)

# Outer constructor for empty object
MolecularDynamicsTrial() = MolecularDynamicsTrial(1,1,0,0,0,0,1,0);

# Comparator for MolecularDynamicsTrial type
function ==( r::MolecularDynamicsTrial, p::MolecularDynamicsTrial )
    for field in names( MolecularDynamicsTrial )
        if r.( field ) != p.( field )
            return false
        end
    end
    return true
end

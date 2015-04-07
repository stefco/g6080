type MolecularState
    # Constants (Mutable, default to Argon/Verlet)
    σ::Real = 3.405e-10             # Length units, from Verlet, in m
    ϵK::Real = 119.8                # Energy units, from Verlet, n Kelvin
    kB::Real = 1.3806e-23           # Boltzmann constant, in J / K
    ϵ::Real = kB * ϵK               # Energy units, from Verlet, in Joules
    mass::Real = 39.948 / 6.0221e23 * 1e-3          # Argon mass, in kg
    τ::Real = sqrt( mass * σ^2 / ( 48ϵ ) )          # Time units, in s
    
    numBodies::Unsigned             # Number of particles
    particleDensity::Real           # Density of particles
    desiredTemperature::Real        # Desired temperature
    dx::Real = 0.032                # Time step; defaults to Verlet's value
    dt::Real = dx * τ               # Time step in seconds
end

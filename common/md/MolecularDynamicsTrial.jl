################################################################################
#
#   Struct for holding all parameters associated with an md simulation
#
################################################################################
abstract MDTrial

################################################################################
#
#   AT SOME POINT, IT WOULD BE GOOD TO HAVE A MACRO FOR AUTOMATICALLY
#   POPULATING MolecularDynamicsTrial TYPE FIELDS
#
################################################################################

################################################################################
#
#   Struct for holding all parameters associated with an md verlet simulation
#
################################################################################
type MolecularDynamicsTrial <: MDTrial
    # Parameters
    numBodies::Int64                # Number of particles
    ρ::Float64                      # Density of particles
    Td::Float64                     # Desired temperature
    dx::Float64                     # Time step; defaults to Verlet's value
    V::Float64                      # Volume of box
    L::Float64                      # Side length of box
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

################################################################################
#
#   Struct for holding all parameters associated with an md metropolis
#   simulation
#
################################################################################
type MDMetropolisTrial <: MDTrial
    # Parameters
    numBodies::Int64                # Number of particles
    ρ::Float64                      # Density of particles
    Td::Float64                     # Desired temperature
    # dx::Float64                     # Time step; defaults to Verlet's value
    V::Float64                      # Volume of box
    L::Float64                      # Side length of box
    steps::Int64                    # How long to run the simulation
    σ::Float64                      # Length scale, from Verlet, in m
    ϵ::Float64                      # Energy scale, in Kelvin
    m::Float64                      # Argon atom mass, in kg
    τ::Float64                      # Time scale, in seconds

    # Dynamic Quantities
    y::Array{Float64, 3}            # Positions
    # v::Array{Float64, 3}            # Velocities
    # f::Array{Float64, 3}            # Accelerations
    pe::Array{Float64, 2}           # Potential Energies
    ke::Array{Float64, 2}           # Kinetic Energies
    pet::Array{Float64, 1}          # Total Potential Energy
    ket::Array{Float64, 1}          # Total Kinetic Energy
    # T::Array{Float64, 1}            # Temperature at each step
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

# Check if trial is finished
function isfinished( r::MolecularDynamicsTrial )
    return r.currentStep == r.steps
end

# Duplicate a trial
function duplicate( r::MolecularDynamicsTrial )
    p = MolecularDynamicsTrial()
    for field in names( MolecularDynamicsTrial )
        p.(field) = r.(field)
    end
    assert(p == r)
    return p
end

# Extract a sub-trial TODO : reimplement using slicedim(A,d,i)
function subtrial( r::MolecularDynamicsTrial, n::Int, m::Int )
    if r.steps < m
        error( "Out of bounds: final step must be less than r.steps" )
    elseif m < n
        error( "first step must be less than final step" )
    elseif n < 1
        error( "Out of bounds: first step must be greater than or equal to 1")
    elseif n > r.currentStep
        error( "first step must be less than current step" )
    else
        steps = m - n + 1
        numBodies = r.numBodies
        p = MolecularDynamicsTrial(numBodies,1,0,0,steps,0,1,0);
        for field in names( MolecularDynamicsTrial )
            if :steps == field
                continue                            # clobber not total steps
            elseif isa( r.(field), Array )
                fsize = size( p.(field) )           # dimensions of the field
                perstep = prod(fsize[1:(end-1)])    # elements per step
                first = (n-1)*perstep + 1           # first element to copy
                last = m*perstep                    # last element to copy
                assert((last-first+1)==perstep*steps)   # correct num elements
                i = 1
                for j in first:last                 # linear indexing to set p
                    p.(field)[i] = r.(field)[j]
                    i+=1
                end
            else
                p.(field) = r.(field)
            end
        end
        p.currentStep = r.currentStep + 1 - n
        return p
    end
end

# Lop off unfinished steps
function mdtruncate( r::MolecularDynamicsTrial )
    if isfinished( r )
        return r
    else
        return subtrial( r, 1, r.currentStep )      # go till last finished step
    end
end

# Add steps, if you need to run for longer
function addsteps( r::MolecularDynamicsTrial, n::Int )
    if n<0
        error( "Must specify non-negative number of steps to add" )
    else
        steps = r.steps + n
        numBodies = r.numBodies
        p = MolecularDynamicsTrial(numBodies,1,0,0,steps,0,1,0);
        for field in names( MolecularDynamicsTrial )
            if :steps == field
                continue                            # clobber not total steps
            elseif isa( r.(field), Array )
                i = 1
                for val in r.(field)[1:end]         # linear indexing to set p
                    p.(field)[i] = val
                    i+=1
                end
            else
                p.(field) = r.(field)               # copy scalar fields
            end
        end
        return p
    end
end

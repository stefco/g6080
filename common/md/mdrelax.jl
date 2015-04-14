################################################################################
#
#   Rescale the velocities of the particles in this trial at the specified step
#
################################################################################
function relax!( r::MolecularDynamicsTrial, n::Int64 )
    r.v[:,:,n] = sqrt(r.Td/r.T[n])*r.v[:,:,n]
end

################################################################################
#
#   Rescale the velocities of the particles at the most recent step
#
################################################################################
function relax!( r::MolecularDynamicsTrial )
    relax!( r, r.currentStep )
end

################################################################################
#
#   Aggressively rescale the velocities to relax the system more quickly
#
################################################################################
function relaxaggresively!( r::MolecularDynamicsTrial, n::Int64, mult::Float64 )
    mult >= 1 || error("multiplier needs must be greater than one.")
    if r.Td >= r.T[n]
        Td = r.Td * mult
    else
        Td = r.Td / mult
    end
    r.v[:,:,n] = sqrt(Td/r.T[n])*r.v[:,:,n]
end

################################################################################
#
#   Aggressively rescale the velocities at most recent step
#
################################################################################
function relaxaggresively!( r::MolecularDynamicsTrial, mult::Float64 )
    relaxaggresively!(r,r.currentStep,mult)
end

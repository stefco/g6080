################################################################################
#
#   Run the velocity Verlet algorithm to calculate steps n through m.
#
#   Rebuild the index before each step.
#   
#   Velocty Verlet works by calculating:
#       1) next position from current velocity and acceleration
#       2) next potential energy and force from next position
#       3) next velocity from current AND next force
#       4) next kinetic energy from next velocity
#
################################################################################
function verlet!(
        r::MolecularDynamicsTrial,
        n::Int64, 
        m::Int64, 
        index::MolecularDynamicsIndex )
    r.currentStep = n-1
    for i in r.currentStep:(m-1)
        index!( r.y[:,:,i], index, r.L )
        r.y[:,:,i+1] = (r.y[:,:,i] + r.v[:,:,i]r.dx + 0.5r.f[:,:,i]r.dx^2) % r.L
        potentialenergyandforce!( r, i+1 )
        r.v[:,:,i+1] = r.v[:,:,i] + 0.5r.f[:,:,i]r.dx + 0.5r.f[:,:,i+1]r.dx
        kineticenergy!( r, i+1 )
        r.e[i+1] = r.ket[i+1] + r.pet[i+1]
        r.currentStep+=1
    end
end

################################################################################
#   Run the Verlet algorithm for another n steps, or until finished
################################################################################
function verlet!( r::MolecularDynamicsTrial, n::Int, i::MolecularDynamicsIndex )
    start = r.currentStep + 1                       # left off on currentStep
    stop = min( r.currentStep + n, r.steps )        # calculate another n steps
    verlet!( r, start, stop, i )
end

################################################################################
#   Run the Verlet algorithm until completion
################################################################################
function verlet!( r::MolecularDynamicsTrial, i::MolecularDynamicsIndex )
    start = r.currentStep + 1
    stop = r.steps
    verlet!( r, start, stop, i )
end


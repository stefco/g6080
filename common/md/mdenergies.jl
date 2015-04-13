################################################################################
#
#   Function for finding the shortest vector twixt two points
#
################################################################################
function Δ( p2::Array{Real,1}, p1::Array{Real,1}, l::Real )
    # Calculate minimum squared distance, accounting for periodicity of the
    # containing box.
    Δp=p2-p1-l
    for i in 1:3
        abs(Δp[i]+l)<abs(Δp[i])&&(Δp[i]+=l)
        abs(Δp[i]+l)<abs(Δp[i])&&(Δp[i]+=l)
    end
    return Δp
end

################################################################################
#
#   Function for finding f⃗, PE, PEtotal, and the Virial at step n
#
#       Calculate kinetic energy by iterating through the interaction list, 
#       calculating rSquared, PE, and F for each interacting body. The 
#       interaction list holds a list of (Int64,Int64) specifying the indices 
#       of the current particle as well as a nearby neighbor.
#   
#       Use Δ(r1,r2,l) to find the shortest distance between pairs, taking 
#       account of periodic boundary conditions. Then calculate dot(Δr⃗,Δr⃗), 
#       which gives the square of the distance between the pairs, which 
#       conveniently is the only value we need for subsequent calculations.
#   
#       Increment f⃗, PE, PEtot, and the Virial by appropriate amounts. Save a 
#       bit of time by holding off on multiplying the virial by 1/2 until all 
#       the contributions to it have been added in.
#
################################################################################
function potentialenergyandforce!( r::MolecularDynamicsTrial, n::Int )
    for pair in interactionlist # TODO : where does interactionlist live?
        Δr⃗ = Δ( r.y[:,pair[1],n], r.y[:,pair[2],n], r.L )
        rsq = dot( Δr⃗, Δr⃗ )
        Δf⃗ = (rsq^-7 - 0.5rsq^-4)Δr⃗
        r.f[:,pair[1],n] += Δf⃗
        r.pe[pair[1],n] += 4(rsq^-6 - rsq^-3)
        r.vir += dot( Δr⃗, Δf⃗ )
    end
    r.pet[n] = sum( r.pe[:,n] )
    r.vir /= 2
end

################################################################################
#
#   Function for finding KE, KEtotal, and T at step n
#
#       In Verlet's units, KE = 24v^2 and T = (2/3)KEavg = (2/3)KEtotal/N.
#
#       1) Calculate KE for each particle
#       2) Sum up the KE to get KEtotal
#       3) Use KEtotal to find the temperature, T
#
################################################################################
function kineticenergy!( r::MolecularDynamicsTrial, n::Int )
    for i=1:r.numBodies
        r.ke[i,n] = 24dot( r.v[:,i,n], r.v[:,i,n] )
    end
    r.ket[n] = sum( r.ke[:,n] )
    r.T[n] = (2/3)r.ket[n] / r.numBodies
end

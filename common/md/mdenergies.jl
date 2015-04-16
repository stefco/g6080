import NumericExtensions

################################################################################
#
#   Type for keeping track of nearby pairs
#
################################################################################
type MolecularDynamicsIndex
    pairs::Array{(Int64,Int64),1}               # tuple vector of nearby pairs
    cells::Array{Array{Int64,1},3}              # grid containers for particles
    δ::Float64                                  # max distance worth including
    λ::Int64                                    # steps before index regen
    Λ::Int64                                    # current step
end

################################################################################
#
#   Convenience constructor for MolecularDynamicsIndex
#
################################################################################
MolecularDynamicsIndex( δ, λ ) =
    MolecularDynamicsIndex( 
        Array((Int64,Int64),0),
        Array(Array{Int64,1},0,0,0),
        δ,
        λ,
        10 )

################################################################################
#
#   Function for finding the shortest vector twixt two points
#   
#       Calculate minimum squared distance, accounting for periodicity of the
#       containing box.
#
################################################################################
function Δ( p2::Array{Float64,1}, p1::Array{Float64,1}, l::Float64, Δp::Array{Float64,1} )
    Δp=p2-p1
    lo2=l÷2
    for i in 1:3
        if Δp[i] > lo2
            Δp[i] -= l
        elseif Δp[i] < lo2
            Δp[i] += l
        end
    end
    return Δp
end

################################################################################
#
#   Simple function for finding nearby pairs
#
#       Iterates through all pairs and, if a pair is sufficiently close, adds
#       it to the list of interacting pairs as a tuple of the particle indices.
#
#       1) If the index counter Λ has reached the planned index lifetime λ,
#          rebuild the index.
#       2) Find the square of the max interaction distance δ; call it δsq
#       3) Clear the list of interacting pairs and reinitialize it as empty.
#       4) Loop through all ordered pairs of particles (j,k), skipping j==k
#       5) Find Δr⃗, the shortest distance between particles j and k
#       6) If the squared magnitude of Δr⃗ is less that δsq, the particles are
#          close enough to interact; push the tuple (j,k) into the pairs list.
#       7) Increment Λ so that it rolls back to 1 on passing λ.
#
################################################################################
function index!( y::Array{Float64, 2}, i::MolecularDynamicsIndex, l::Float64 )
    if i.Λ == i.λ
        print("Rebuilding the index... ")
        δsq = i.δ^2
        i.pairs = Array((Int64,Int64),0)
        for j in 1:size(y)[2]
            for k in 1:size(y)[2]
                if j == k
                    continue
                else
                    Δr⃗ = Δ(y[:,j],y[:,k],l)
                    if δsq >= dot(Δr⃗,Δr⃗)
                        push!( i.pairs, (j,k) )
        end; end; end; end
        print("done. ")
    end
    i.Λ = ( i.Λ % i.λ ) + 1
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
function potentialenergyandforce!( r::MolecularDynamicsTrial, n::Int64 )
    m = r.numBodies
    l = r.L
    x = r.y[:,:,n]

    ΔR⃗ = reshape(x,3,m) .- reshape(x,3,1,m) # All interatomic distances, aperiodic
    ΔR⃗ *= (1/l)                             # Scale out box length and shift to
    ΔR⃗ += 0.5                               # middle, so that 0 < ΔR⃗min < 1.
    ΔR⃗ = mod(ΔR⃗,1)                          # Pick values that match this criterion.
    ΔR⃗ -= 0.5                               # Shift values back into place, and
    ΔR⃗ *= l                                 # scale them back up.

    Rn2 = NumericExtensions.sumsq(ΔR⃗,1)     # All the minimum square distances;
    Rn2 = 1 ./ Rn2                          # We will use Rn2 = R^-2, hence the name
    for i in 1:m                            # Avoid infinite self-interactions
        Rn2[1,i,i] = 0.0
    end
    Rn6 = Rn2 .* Rn2 .* Rn2                 # Rn6 = Rn^-6, which will be handy

    ΔPE = 4(Rn6 .* Rn6 - Rn6)
    Δvir = 0.25Rn6 + 0.125ΔPE
    Δf⃗ = 2Rn2 .* (Δvir .* ΔR⃗)

    r.pe[:,n] = sum(ΔPE,2)
    r.pet[n] = sum(r.pe[:,n])
    r.vir[n] = sum(Δvir)
    r.f[:,:,n] = sum(Δf⃗,3)
end

################################################################################
#
#   Calculate just PE and the virial
#
################################################################################
function potentialenergy!( r::MolecularDynamicsTrial, n::Int64 )
    m = r.numBodies
    x = r.y[:,:,n]

    x2 = NumericExtensions.sumsq(x,1)       # length of each column in x
    Rn2 = x2 .+ reshape(x2,m,1)             # xi^2 +xj^2
    gemm!('T', 'N', -2.0, x, x, 1.0, Rn2)   # add -2xi'xj to Rn2

    Rn2 = 1 ./ Rn2                          # Rn2 = R^-2, hence the name
    for i in 1:m                            # No self-interaction
        Rn2[i,i] = 0.0
    end
    Rn6 = Rn2 .* Rn2 .* Rn2                 # Rn6 = R^-6

    ΔPE = 4(Rn6 .* Rn6 - Rn6)                # Multiply by 4 after getting virial
    Δvir = 0.25Rn6 + 0.125ΔPE                     # Virial contribution

    r.pe[:,n] = sum(ΔPE,1)
    r.pet[n] = sum(r.pe[:,n])
    r.vir[n] = sum(Δvir)
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
function kineticenergy!( r::MolecularDynamicsTrial, n::Int64 )
    print("Calculating kinetic energies... ")
    for i=1:r.numBodies
        r.ke[i,n] = 24dot( r.v[:,i,n], r.v[:,i,n] )
    end
    r.ket[n] = sum( r.ke[:,n] )
    r.T[n] = (2/3)r.ket[n] / r.numBodies
    print("done. ")
end

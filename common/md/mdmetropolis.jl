import NumericExtensions.sumsq

################################################################################
#
#   METROPOLIS can take a full MolecularDynamicsTrial object, or the reduced
#   MDMetropolisTrial object, which is far more space-efficient, since the
#   fields of the latter are a subset of the former.
#
#   Otherwise works much like verlet!(r,n,m,δ), i.e. calculates values for steps
#   n through m on the MDTrial r.
#
#   The Metropolis algorithm itself works as follows:
#       1) Loop through particle positions, creating a series of proposed 
#          random moves a distance δ from previous positions
#       2) calculate new particle energy at this position
#       3) calculate 
#
################################################################################
function metropolis!( r::MDTrial, n::Int64, m::Int64, δ::Float64 )
    # Parameters
    N = r.numBodies
    l = r.L
    oneol = 1/l                             # Use this so much; set it now
    β = 1 / (r.Td)                          # T, E in units of ϵ, so kB cancels
    r.currentStep = n-1

    # Allocations
    A = Array(Float64, 3, N)                # all accepted new positions
    P = Array(Float64, 3, N)                # all proposed positions
    
    p⃗ = Array(Float64, 3)                   # a proposed position
    Δp⃗ = Array(Float64, 3, N)               # dist btwn p⃗ and all Y's

    pPE = 0.0                               # proposed PE for a particle
    ΔPE = 0.0                               # ΔPE resulting from a proposed step
    
    α = 0.0                                 # prob of switching to higher pPE
    acc = 0.0                               # accepted movements per step
    
    R2 = Array(Float64, 1, N)               # sq distances betwn Δp⃗ and Y's
    Rn2 = Array(Float64, 1, N)              # 1./R2
    Rn6 = Array(Float64, 1, N)              # Rn2^3
    
    # Calculations
    for i in r.currentStep:(m-1)
        acc = 0.0
        print("START Metropolis step ",r.currentStep+1,". ")
        A = r.y[:,:,i]
        P = A - (2.0rand(3,N) - 1.0)δ
        print("Running acc/rej... ")
        for j in 1:N
            # Handle periodic boundaries, find shortest interparticle Δp⃗
            p⃗ = P[:,j]
            Δp⃗ = p⃗ .- A
            Δp⃗ *= oneol
            Δp⃗ += 0.5
            Δp⃗ = mod(Δp⃗,1)
            Δp⃗ -= 0.5
            Δp⃗ *= l

            # Find squared distances and compute energy
            R2 = sumsq(Δp⃗,1)
            Rn2 = 1 ./ R2
            Rn2[j] = 0                      # avoid self-interaction
            Rn6 = Rn2 .* Rn2 .* Rn2

            # Find ΔPE due to the proposed move
            pPE = 4.0(sumsq(Rn6,2) - sum(Rn6))
            ΔPE = 2.0(pPE - r.pe[j,i])

            # If the new state has lower energy, accept new position
            @assert length(ΔPE) == 1
            if ΔPE[1] < 0.0
                A[:,j] = p⃗                  # accept new position
                acc += 1
            else
                α = exp( -β*ΔPE[1] )
                if α > rand()
                    A[:,j] = p⃗              # accept new position
                    acc += 1
                end
            end
        end
        r.y[:,:,i+1] = A                    # accepted values become new r.y
        print("done with acc=",acc/N,". Find PE and virial... ")
        potentialenergy!(r, i+1)            # get pe, pet, and virial
        r.currentStep += 1
        println("done.")
    end
end

function metropolis!(r::MDTrial, n::Int64, δ::Float64)
    start = r.currentStep + 1
    stop = min( r.currentStep + n, r.steps )
    metropolis!( r, start, stop, δ )
end

function metropolis!(r::MDTrial, δ::Float64)
    start = r.currentStep + 1
    stop = r.steps
    metropolis!( r, start, stop, δ )
end

# Wolff single-cluster flip code
# Stefan Countryman
# stc2117@columbia.edu

type SpinArray
    Nx::Int64
    Ny::Int64
    N ::Int64
    σ ::Array{Int64,2}
end

# Initialize a spin-up spin array with specified dimensions
SpinArray(Nx::Int64, Ny::Int64) = SpinArray(Nx, Ny, Nx*Ny, ones(Int64, Ny, Nx))

function magnetization(A::SpinArray)
    return sum(A.σ)/A.N
end

function flip!(x::Int64, y::Int64, A::SpinArray, Δσ::Array{Int64,2}, flips::Int64, p::Float64)
    
    s0 = Δσ[x,y]                                # current spin
    Δσ[x,y] = -s0                               # flip this spot
    flips += 1                                  # add a flip
    X = [mod1((x+1), A.Nx), mod1((x-1), A.Nx)]  # x-coordinates of neighbors
    Y = [mod1((y+1), A.Ny), mod1((y-1), A.Ny)]  # y-coordinates of neighbors
    for xx in X
        for yy in Y
            if (s0 == Δσ[xx,yy] && rand() < p)
                flips = flip!(xx, yy, A, Δσ, flips, p)
            end
        end
    end
    return flips
end

function cluster(Nx::Int64, Ny::Int64, steps::Int64, J::Float64, B::Float64, T::Float64)

    print("Finding normalized magnetization for J=$J, B=$B, T=$T... ")
    p = 1-exp(-2J/T)
    A = SpinArray(Nx, Ny)
    Δσ = zeros(A.σ)                             # proposed spins
    magnetizations = zeros(steps)
    attempts = 0
    
    for n in 1:steps

        stilltrying = true
        while stilltrying
            # another attempt
            attempts += 1

            # randomly pick a spot to flip
            x = int(ceil(rand() * Nx))
            y = int(ceil(rand() * Ny))
    
            # reset variables
            Δσ[:] = A.σ[:]
            s0 = Δσ[x,y]
            flips = 0
    
            # MUST RECALCULATE ALL ENERGIES
            # propose flip
            flips = flip!(x, y, A, Δσ, flips, p)
            
            # accept/reject
            if rand() < exp(-2B*s0*flips/T) 
                (A.σ[:] = Δσ[:])
                stilltrying = false
            end
        end

        # get magnetization
        magnetizations[n] = magnetization(A)
        # n%100 == 0 && println("\tMagnetization for step $n: ",magnetizations[n])
    end
    
    println("done with ",(steps/attempts)," accept ratio.")
    return magnetizations
end

function M(T,J)
    return (1 - (sinh(2J./T)).^-4).^(1/8)
end



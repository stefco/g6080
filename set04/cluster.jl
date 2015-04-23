# Wolff single-cluster flip code
# Stefan Countryman
# stc2117@columbia.edu

type SpinArray
    Nx::Int64
    Ny::Int64
    N ::Int64
    σ ::Array{Int64,2}
    Eold::Float64
end

# Initialize a spin-up spin array with specified dimensions
SpinArray(Nx::Int64, Ny::Int64) = SpinArray(Nx, Ny, Nx*Ny, ones(Int64, Ny, Nx), 0.0)

function magnetization(A::SpinArray)
    return sum(A.σ)/A.N
end

function flip!(x::Int64, y::Int64, A::SpinArray, Δσ::Array{Int64,2}, flips::Int64, p::Float64)
    
    s0 = Δσ[x,y]                                # current spin
    Δσ[x,y] = -s0                               # flip this spot
    flips += 1                                  # add a flip
    nei = [
        (x, mod1(y+1, A.Ny))
        (x, mod1(y-1, A.Ny))
        (mod1(x+1, A.Nx), y)
        (mod1(x-1, A.Nx), y)
    ]

    for (xx, yy) in nei
        if (s0 == Δσ[xx,yy] && rand() < p)
            flips = flip!(xx, yy, A, Δσ, flips, p)
        end
    end
    return flips
end

# energies are -2J*∑s1*s2/T; sum them all up
function interactionenergy(spins::Array{Int64,2}, J::Float64)
    interactions = 0
    horshift = [spins[:,2:end]  spins[:,1]]
    horshift .*= spins
    vershift = [spins[2:end,:]; spins[1,:]]
    vershift .*= spins
    interactions += sum(horshift)
    interactions += sum(vershift)
    interactions *= -2
    return J * float64(interactions)
end

function magneticenergy(spins::Array{Int64,2}, B::Float64)
    spin = sum(spins)
    spin *= -1
    return B * float(spin)
end

function cluster(Nx::Int64, Ny::Int64, steps::Int64, J::Float64, B::Float64, T::Float64)

    print("Finding normalized magnetization for J=$J, B=$B, T=$T... ")
    p = 0.46 # 1-exp(-2J/T)
    A = SpinArray(Nx, Ny)
    A.Eold = interactionenergy(A.σ, J) + magneticenergy(A.σ, B)
    Δσ = zeros(A.σ)                             # proposed spins
    magnetizations = zeros(steps)
    accepts = 0
    Enew = 0.0
    ΔE = 0.0
    Pacc = 0.0
    
    for n in 1:steps

        # randomly pick a spot to flip
        x = int(ceil(rand() * Nx))
        y = int(ceil(rand() * Ny))
    
        # reset variables
        Δσ[:] = A.σ[:]
        s0 = Δσ[x,y]
        flips = 0
    
        # propose flip
        flips = flip!(x, y, A, Δσ, flips, p)
        showspins(Δσ)

        # calculate energies; ΔE = (Eint + Emag) - Eold
        Enew = interactionenergy(Δσ, J)
        Enew += magneticenergy(Δσ, B)
        ΔE = Enew - A.Eold
        Pacc = exp(-ΔE/T)  
            
        # accept/reject
        if rand() < Pacc
            A.σ[:] = Δσ[:]
            A.Eold = Enew
            println("Accepted step $n")
            accepts += 1
        end

        # get magnetization
        magnetizations[n] = magnetization(A)
        # n%100 == 0 && println("\tMagnetization for step $n: ",magnetizations[n])
    end
    
    println("done with ",(accepts/steps)," accept ratio.")
    return magnetizations, accepts/steps
end

function M(T,J)
    return (1 - (sinh(2J./T)).^-4).^(1/8)
end

function showspins(spins::Array{Int64,2})
    up = "██"
    down = "  "
    println("Spins (\"$up\" is ↑, \"$down\" is ↓):")
    for nrow in 1:size(spins)[2]
        for spin in spins[:,nrow]
            if spin == 1
                print(up)
            elseif spin == -1
                print(down)
            else
                error("Spin array must have values of 1 or -1.")
            end
        end
        println()
    end
end
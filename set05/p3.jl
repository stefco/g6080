module p3

# Return boundary conditions values and indices, ϕb and β, for problem 2, 
# with N subdivisions
function boundaryconditions3(N::Integer)
    ϕb = zeros(N+1,N+1)                 # boundary conditions
    β  = zeros(Int,N+1,N+1)             # boundary condition locations

    # Conditions on the perimeter
    ϕb[:,1] = ϕb[:,N+1] = ϕb[1,:] = ϕb[N+1,:] = 0.0
    β[:,1]  = β[:,N+1]  = β[1,:]  = β[N+1,:]  = 1

    # Conditions on the interior
    # xmin < x < xmax  →  N*xmin + 1 < i < N*xmax + 1
    xmin = 0.5
    xmax = 0.75
    ymin = 0.625
    ymax = 0.875

    xs = int([ceil(N*xmin+1):floor(N*xmax+1)])
    ys = int([ceil(N*ymin+1):floor(N*ymax+1)])

    ϕb[xs,ys] = 0.0
    β[xs,ys]  = 1

    return ϕb, β
end

# Calculate the potential
function V(x::Real, y::Real, N::Integer, q::Real)

    a = b = 1:(N+1)
    a .-= (N*x + 0.5)           # avoid singularity
    b .-= (N*y + 0.5)
    a .*= a
    b .*= b
    return -N*q ./ √(a .+ b')

end

V(N::Integer, q::Real) = V(0.5, 0.25, N, q)

function minor(A::Matrix, β::Vector)
    
    size(A)[1] == size(A)[2] || error("Matrix must be square")

    N = sum(1 .- β)             # dimensions of new matrix

    B = Array(typeof(A[1]),N,N)

    ii = 1
    jj = 1

    for j in 1:size(A)[1]
        β[j] == 1 && continue
        for i in 1:size(A)[2]
            β[i] == 1 && continue
            B[ii,jj] = A[i,j]
            ii += 1
        end
        jj +=1
        ii = 1
    end

    return B
end

# repopulate a vector v with zeros for fixed parameters
function addwater(v::Vector, β::Vector)

    vf = Array(typeof(v[1]),length(β))
    ii = 1

    length(v) == sum(1-β) || error("v and β have different num free params")

    for i in 1:length(β)
        if β[i] == 1
            vf[i] = 0               # put zero for a fixed value
        elseif β[i] == 0
            vf[i] = v[ii]           # fill in the free value
            ii += 1                 # prep the next free value
        else
            error("β can only have 0 or 1 values")
        end
    end
    return vf

end

end

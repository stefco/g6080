module p2

# Return boundary conditions values and indices, ϕb and β, for problem 2, 
# with N subdivisions
function boundaryconditions2(N::Integer)
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

# Create the inverse iteration matrix, A = Δ - μI, of size N
function invit(Δ::Array, μ::Real, β::Array)
    return -5Δ - μ .* eye(size(Δ)[1]) .* ((1.-β)*(1.-β)')
end

end

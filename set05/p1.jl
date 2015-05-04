module p1

# Return boundary conditions values and indices, ϕb and β, for problem 1, 
# with N subdivisions
function boundaryconditions1(N::Integer)
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

    ϕb[xs,ys] = 10.0
    β[xs,ys]  = 1

    return ϕb, β
end

# Return the charge distribution for problem 1
function chargedensity(N::Integer, x::Real, y::Real, q::Real)
    ρ = zeros(N+1,N+1)

    x <= 1.0 || error("x must be less than 1")
    y <= 1.0 || error("y must be less than 1")

    # Spread the point charge across four points; treat it like uniform density
    ρ[floor(x*N+1),floor(y*N+1)] = N*N*q/4
    ρ[floor(x*N+2),floor(y*N+1)] = N*N*q/4
    ρ[floor(x*N+1),floor(y*N+2)] = N*N*q/4
    ρ[floor(x*N+2),floor(y*N+2)] = N*N*q/4

    return ρ
end

chargedensity1(N::Integer) = chargedensity(N, 0.5, 0.25, 20)

# Get the grid location nearest to a given x,y value
function nearestpoint(x::Real, y::Real, N::Integer)
    
    i = int(x*N+1)
    j = int(y*N+1)

    return i, j

end

end

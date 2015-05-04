module Laplacian

# Create a Kronecker Delta matrix with the given relationship and dimensions
function δ(fi::Function, fj::Function, n::Integer)
    r = zeros(n,n)
    for j in [1:n]
        for i in [1:n]
            fi(i) == fj(j) && (r[i,j] = 1.0)
        end
    end
    return r
end

kronecker = δ

# Create a laplacian matrix with specified grid size and grid divisions
function lap(N::Integer, L::Real)

    Np1sq = (N+1)*(N+1)                 # have (N+1)^2 grid points
    a = -4.0δ(i->i, j->j, Np1sq)
    a += δ(i->mod1(i+1, Np1sq), j->j, Np1sq)
    a += δ(i->mod1(i-1, Np1sq), j->j, Np1sq)
    a += δ(i->mod1(i+N+1, Np1sq), j->j, Np1sq)
    a += δ(i->mod1(i-N-1, Np1sq), j->j, Np1sq)
    return a * (N*N / (L*L))

end

lap(N::Integer) = lap(N,1)

# Split a laplacian matrix into free and boundary parts
function opsplit(A::Array, β::Array{Int,1})

    ndims(A) == 2 || error("operator matrix must be 2-d")
    length(β) == size(A)[1] || error("specify which indices are boundary conditions")
    size(A)[1] == size(A)[2] || error("operator matrix must be square")

    B = zeros(size(A))
    F = zeros(size(A))
    # Set the boundary condition columns equal to the corresponding columns of A
    for i in 1:length(β)
        β[i] == 1 ? B[:,i] = A[:,i] : F[:,i] = A[:,i]
    end
    return F, B

end

# β[i] = 1 iff ϕ[i] is a boundary condition
function rhoeff(ρ::Array, Δb::Array, ϕb::Array, β::Array{Int,1})
    return (ρ + Δb*ϕb) .* (1 .- β)
end

# Cancel out bound rows to get Δeffective
function deltaeff(Δf, β::Array{Int,1})
    return Δf .* (1 .- β)
end

end

module Laplacian

# Cyclically permute slices of the array in direction of increasing index
function cycle(A::Array, n::Int64, dim::Int64)

    l = size(A)[dim]
    m = mod1(n, l)
    return cat(dim, slicedim(A, dim, (m+1):l), slicedim(A, dim, 1:m))

end

# Compute the discrete Laplacian at every point in a grid
function ∇²(F::Array{Float64, 2}, h::Float64)
    
    A = 4F
    A -= cycle(A,  1, 1)        # neighbor on right
    A -= cycle(A, -1, 1)        # neighbor on left
    A -= cycle(A,  1, 2)        # neighbor above
    A -= cycle(A, -1, 2)        # neighbor below

    c = 1 / (h * h)             # constant for grid spacing
    return A .* c                # multiplication is faster than division

end

# Create a Kronecker Delta matrix with the given relationship and dimensions
function δ(fi::Function, fj::Function, n::Integer)
    r = zeros(Int64,n,n)
    for j in [1:n]
        for i in [1:n]
            fi(i) == fj(j) && (r[i,j] = 1)
        end
    end
    return r
end

# Create a laplacian matrix with specified grid size and grid divisions
function lap(N::Integer, L::Real)

    Nsq = N*N                       # have N^2 grid points
    a = 4δ(i->i, j->j, Nsq)
    a -= δ(i->mod1(i+1, Nsq), j->j, Nsq)
    a -= δ(i->mod1(i-1, Nsq), j->j, Nsq)
    a -= δ(i->mod1(i+N, Nsq), j->j, Nsq)
    a -= δ(i->mod1(i-N, Nsq), j->j, Nsq)
    return a * (Nsq / (L*L))

end

lap(N::Integer) = lap(N,1)

end

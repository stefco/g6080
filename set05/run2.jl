import Laplacian
import ConjugateGradient
import p2

# Run question 1 with grid size N
function q2(N::Integer, μs::Array)

    # Boundary conditions as matrices
    Ψb, B = p2.boundaryconditions2(N)

    # Reshape conditions into vectors
    ψb = reshape(Ψb, (N+1)^2)
    β  = reshape(B,  (N+1)^2)

    # Laplacian matrix
    Δ = Laplacian.lap(N)
    Δf, Δb = Laplacian.opsplit(Δ, β)
    Δe = Laplacian.deltaeff(Δf, β)

    # First guess should be zero for boundary points
    b0 = ones(Int64,(N+1)^2) .* (1 .- β)  / (N+1 - sum(β))

    # Compute Bμ (in parallel, of course)
    Bs = Array(Float64, length(μs))
    @sync @parallel for i in 1:length(μs)
        # Inverse iteration matrix A
        println("On μ = $(μs[i])")
        A = p2.invit(Δe, μs[i], β)

        # Grip it and rip it
        Bs[i] = norm(ConjugateGradient.conjgrad(A, b0))
    end

    # return the B values
    return Bs, Δe

end

# Return X, Y, and Z coordinates of ϕ gridi matrix, assuming 1×1 grid size
function xyzcoords(ϕ)

    N = size(ϕ)[1] - 1
    a = 1/N

    ii = [0:N]  .* ones(N+1)'
    jj = [0:N]' .* ones(N+1)

    xs = a .* ii
    ys = a .* jj
    zs = ϕ

    X = reshape(xs, (N+1)^2)
    Y = reshape(ys, (N+1)^2)
    Z = reshape(zs, (N+1)^2)
    
    return X, Y, Z

end

function writexyxcoords(prefix,ϕ)

    X, Y, Z = xyzcoords(ϕ)

    writedlm("$prefix-phi-X.dat", X)
    writedlm("$prefix-phi-Y.dat", Y)
    writedlm("$prefix-phi-Z.dat", Z)

end

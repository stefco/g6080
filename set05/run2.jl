import Laplacian
import ConjugateGradient
import p2

# Run question 2 with grid size N
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
    b0 = ones(Int64,(N+1)^2) .* (1 .- β)  / √sum(1.-β)

    # Compute Bμ (in parallel, of course)
    Bs = Array(Float64, length(μs))
    @sync @parallel for i in 1:length(μs)
        # Inverse iteration matrix A
        println("On μ = $(μs[i])")
        A = p2.invit(-5Δe, μs[i], β)

        # Grip it and rip it
        b = ConjugateGradient.conjgrad(A, b0)
        Bs[i] = dot(b,b)
    end

    # return the B values
    return Bs, Δe

end

# Find eigenvectors corresponding to our eigenvalues
function q2vectors(N::Integer, μ::Float64)

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
    ψ0 = ones(Int64,(N+1)^2) .* (1 .- β)  / √sum(1.-β)

    # Grip it and rip it; calculate ϕfree, use it to find total ψ
    A = p2.invit(-5Δe, μ, β)
    ψf = ConjugateGradient.conjgrad(A, ψ0)
    ψf /= norm(ψf)
    ψ = β.*ψb + (1.-β).*ψf

    # Reshape ψ to give the physical grid and return that value
    return reshape(ψ, N+1, N+1)

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

import Laplacian
import ConjugateGradient
import p1

# Run question 1 with grid size N
function q1(N::Integer)

    # Boundary conditions as matrices
    Vb, B = p1.boundaryconditions1(N)
    @show r = p1.chargedensity1(N)

    # Reshape conditions into vectors
    ϕb = reshape(Vb, (N+1)^2)
    β  = reshape(B,  (N+1)^2)
    ρ  = reshape(r,  (N+1)^2)

    # Laplacian matrix
    Δ = Laplacian.lap(N)
    Δf, Δb = Laplacian.opsplit(Δ, β)
    Δe = Laplacian.deltaeff(Δf, β)

    # Effective charge density
    @show ρe = Laplacian.rhoeff(ρ, Δb, ϕb, β)

    # Grip it and rip it; calculate ϕfree, use it to find total ϕ
    ϕf = ConjugateGradient.conjgrad(Δe, -ρe)
    ϕ = β.*ϕb + (1.-β).*ϕf

    # Reshape ϕ to give the physical grid and return that value
    return reshape(ϕ, N+1, N+1), r, reshape(ρe, N+1, N+1)

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

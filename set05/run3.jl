import Laplacian
import p3

# Generate a clean Hamiltonian
function q3ham(N::Integer, q::Real)

    # Boundary conditions as matrices
    Ψb, B = p3.boundaryconditions3(N)

    # Reshape conditions into vectors
    ψb = reshape(Ψb, (N+1)^2)
    β  = reshape(B,  (N+1)^2)
    V  = reshape(p3.V(N,q), (N+1)^2)

    # Hamiltonian matrix
    H = -5Laplacian.lap(N) - V.*eye((N+1)^2)
    Hf, Hb = Laplacian.opsplit(H, β)
    He = Laplacian.deltaeff(Hf, β)
    Hm = p3.minor(He, β)
    return Hm, He, ψb, β

end

function q3groundeig(Hm::Matrix, steps::Integer, tolerance::Real)

    v = rand(size(Hm)[1])
    v /= norm(v)

    # using the minor matrix kill most zero-valued eigenvalues
    T, K = Lanczos.lanczos(Hm,v,steps)

    # find smallest (nonzero) eigenvector
    D, V = eig(T)
    k = findnext(x -> >(x,1.0), D, 1)


    # old eigenvalue
    λ = D[k]
    err = 100.0

    # best guess at small eigenvalue
    vold = v
    v = K*V[:,k]
    v /= norm(v)

    while err > tolerance
        T, K = Lanczos.lanczos(Hm,v,steps)
        D, V = eig(T)
        k = findnext(x -> >(x,1.0), D, 1)

        # how far off our eigenvalue guess is from the real thing
        err = norm((Hm - λ*eye(size(Hm)[1]))*v)
        λ = D[k]                    # smallest nonzero eigenvalue

        vold = v
        v = K*V[:,k]
        v /= norm(v)

    end

    # return the eigenvalue and eigenvector
    return λ, v

end

# Get the ground state energy and wavefunction, with the latter as a matrix
function q3groundpsi(N::Integer, q::Real, steps::Integer, tolerance::Real)
    Hm, He, ψb, β = q3ham(N, q)
    E0, ψ0 = q3groundeig(Hm, steps, tolerance)
    return E0, reshape(p3.addwater(ψ0,β) + ψb, N+1, N+1)
end

q3groundpsi(N::Integer, q::Real) = q3groundpsi(N, q, 10, 1e-3)

# Return ground state energies and wavefunctions for different q values
function q3groundpsi(N::Integer, qs::Vector)

    Es = Array(Float64, length(qs))
    ψs = Array(Float64, N+1, N+1, length(qs))

    @sync @parallel for i in 1:length(qs)
        Es[i], ψs[:,:,i] = q3groundpsi(N,qs[i])
    end

    return Es, ψs

end

# Find probability that we're in the lower half
function q3P(N::Integer, q::Real, steps::Integer, tolerance::Real)
    E0, ψ0 = q3groundpsi(N, q, steps, tolerance)
    r = ψ0[:,1:int(floor(N/2+1))]   # just bottom half
    return sum(r.*r)                # ψ*ψ
end

# Something we can feed to an equation solver
q3P(q::Real) = q3P(32, q, 20, 1e-3)

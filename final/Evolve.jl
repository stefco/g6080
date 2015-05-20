module Evolve

import Laplacian

# Return the x values associated with a given choice of grid divisions and grid
# width
function X(L :: Real, N :: Integer, xmin :: Real)
    return (L/N)*[0:N]+xmin
end
X(L :: Real, N :: Integer) = X(L, N, 0)

# Return a time evolution matrix A for a given potential energy vector and box width.
function A(V :: Vector, L :: Real)
    
    ϵ = L / (length(V) - 1)

    AA  = SymTridiagonal(0.25(1+ϵ^2*V)im+0.5,-0.125im*ones(length(V)-1))
    AA.dv[1] = AA.dv[end] = 1
    AA.ev[1] = AA.ev[end] = zero(typeof(AA[1,1]))
    return AA

end
A(V :: Function, L :: Real, N :: Integer, xmin :: Real) = A(V(X(L,N,xmin)),L)
A(V :: Function, L :: Real, N :: Integer) = A(V, L, N, 0)

# Return a reverse time evolution matrix B for a given potential energy vector 
# and box width.
function B(V :: Vector, L :: Real)
    
    ϵ = L / (length(V) - 1)

    AA  = SymTridiagonal(0.25(1+ϵ^2*V)im-0.5,-0.125im*ones(length(V)-1))
    AA.dv[1] = AA.dv[end] = 1
    AA.ev[1] = AA.ev[end] = zero(typeof(AA[1,1]))
    return AA

end
B(V :: Function, L :: Real, N :: Integer, xmin :: Real) = B(V(X(L,N,xmin)),L)
B(V :: Function, L :: Real, N :: Integer) = B(V, L, N, 0)

# Return a constructor for an update matrix A with the specified potential energy
# function
function A(V :: Function)
    return (L :: Real, N :: Integer, xmin :: Real) -> A(V, L, N, xmin)
end

# Return the Hamiltonian for a given potential energy vector and box width.
function H(V :: Vector, L :: Real)

    ϵ = L / (length(V) - 1)

    AA = SymTridiagonal((ϵ^-2)-V, -0.5(ϵ^-2)ones(length(V)-1))
    AA.dv[1] = AA.dv[end] = 1
    AA.ev[1] = AA.ev[end] = zero(typeof(AA[1,1]))
    return AA

end
H(V :: Function, L :: Real, N :: Integer, xmin :: Real) = H(V(X(L,N,xmin)),L)
H(V :: Function, L :: Real, N :: Integer) = H(V, L, N, 0)

# Method for setting ϕ₀
function Φ(k₀ :: Real, σ :: Real, x₀ :: Real)

    return (L :: Real, N :: Integer, xmin :: Real) -> begin
        xs = (L/N) .* [0:N] + xmin
        swp = map((x) -> exp(im*k₀*x - (x-x₀)^2 / (2σ^2)), xs)
        swp *= (1 / (π*σ^2))^(1/4)          # normalize
        swp[1] = swp[end] = 0               # enforce ∞ square well boundary
        return swp
    end

end
Φ(k₀ :: Real, σ :: Real) = Φ(k₀, σ, 0)          # start in the center

# Function for evolving a given state vector through n additional steps 
# (first step not in result). Keep running until done or until test fails.
function evolve(A :: AbstractArray, ϕ :: Vector, n :: Integer, test :: Function)

    ϕs = zeros(Complex, length(ϕ), n)

    # Solve for χ using Thomas algorithm
    χ  = tommy(A, ϕ)
    ϕs[:,1] = χ - ϕ

    for j in 1:(n-1)
        if test(ϕs[:,j]) 
            ϕs = ϕs[:,1:j]
            break
        end
        χ = tommy(A, ϕs[:,j])
        ϕs[:,j+1] = χ - ϕs[:,j]
    end

    return ϕs

end
evolve(A :: AbstractArray, ϕ :: Vector, n :: Integer) = evolve(A, ϕ, n, x->true)

# Method for evolving a state vector by passing functions
evolve(A :: Function, ϕ :: Function, n :: Integer, L :: Real, N :: Integer, xmin :: Real) = evolve(A(L,N,xmin), ϕ(L,N,xmin), n)

# Solve Ax = b where A is tridiagonal using the Thomas algorithm
# based on MATLAB code: http://www.mathworks.com/matlabcentral/fileexchange/1359-fast-tridiagonal-system-solver/content/thomas.m
function tommy(a :: Vector, b :: Vector, c :: Vector, d :: Vector)

    n = length(a)
    m = zeros(typeof(a[1]), n)
    l = zeros(typeof(a[1]), n-1)
    y = zeros(typeof(a[1]), n)
    x = zeros(typeof(a[1]), n)

    m[1] = a[1]
    y[1] = d[1]

    # Perform LU factorization and forward substitution to solve Ly = d for y
    for i in 2:n
        l[i-1] = c[i-1] / m[i-1]
        m[i] = a[i] - l[i-1] * b[i-1]

        y[i] = d[i] - l[i-1] * y[i-1]
    end

    x[n] = y[n] / m[n]

    # Begin backward substitution to solve Ux = y for x
    for i in (n-1):-1:1
        x[i] = (y[i] - b[i] * x[i+1]) / m[i]
    end

    return x

end

tommy(A :: Matrix, d :: Vector) = tommy(diag(A), diag(A, -1), diag(A, +1), d)
tommy(A :: Tridiagonal, d :: Vector) = tommy(A.d, A.du, A.dl, d)
tommy(A :: SymTridiagonal, d :: Vector) = tommy(A.dv, A.ev, A.ev, d)

# Create a pocket function centered at zero
function pocket(L :: Real, h :: Real , b :: Real, s :: Real, w :: Real, xmin :: Real, x :: Real)
    if x < (xmin + 0.5L - 0.5w - s) || x > (xmin + 0.5L + 0.5w +s)
        V = 0
    elseif x < (xmin + 0.5L - 0.5w) || x > (xmin + 0.5L + 0.5w)
        V = h
    else
        V = b
    end
    return V
end
 
end

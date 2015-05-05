module Lanczos

function lanczos(A,v1,steps)

    size(A)[1] == size(A)[2] || error("matrix must be square")

    ∘ = dot

    v = zeros(size(v1))
    r = v1
    b = 1

    T = Array(Float64,steps+1,steps+1)      # Lanczos matrix
    K = Array(Float64,size(A)[1],steps+1)   # vectors in Krylov space

    # First round
    vold = v
    v = r/b
    Av = A*v
    a = v∘Av
    T[1,1] = a
    r = Av - v*a - vold*b

    # Loop
    for i in 1:steps
        b = norm(r)
        T[i,i+1] = T[i+1,i] = b
        vold = v
        v = r/b                     # next v
        Av = A*v
        a = v∘Av
        T[i+1,i+1] = a
        r = Av - v*a - vold*b       # orthogonalize
    end

    return T, K

end

lanczos(A,v1) = lanczos(A,v1,size(A)[1])

end

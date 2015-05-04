module ConjugateGradient

function conjgrad(A::Array, b::Array, N::Integer, x1::Array)

    issym(A) || error("Must provide symmetric matrix")

    ∘ = dot                                 # convenience notation

    x = x1
    p = r = b - A*x
    rrold = r∘r

    for i in 1:N
        Ap = A*p
        α = rrold / (p∘Ap)                  # scale correction
        x += α*p                            # add correction
        r -= α*Ap
        rr = r∘r
        √rr < 1e-16 && break                # stop if answer is close
        p *= (rr / rrold)
        p += r
        rrold = rr
    end 

    return x

end

# convenience methods
conjgrad(A, b, N::Integer) = conjgrad(A, b, N, ones(size(b)))
conjgrad(A, b, x1::Array) =  conjgrad(A, b, length(b), x1)
conjgrad(A, b) = conjgrad(A, b, length(b), ones(size(b)))

end

function autocorrelation(M, v, a, ns)
    correlations = Float64[];
    # we can take vector n values as input
    for n in ns
        # Sum over M-n terms (v[i+n,a]-v̂[a])(v[i,a]-v̂[a]), divide by (M-n),
        #   is the same as taking the mean of those terms for i in [1:M-n].
        push!(correlations, mean(
            ( v[(n+1):M,a] - v̄̂[a] ) .* ( v[1:(M-n),a] - v̄̂ [a] )
        ));
    end
    return correlations;
end

function autocorrelation(M, v, a, ns)
    correlations = Float64[];
    # we can take vector n values as input
    for n in ns
        push!(correlations, mean(
            ( v[(n+1):M,a] - v̄̂[a] ) .* ( v[1:(M-n),a] - v̄̂ [a] )
        ));
    end
    return correlations;
end

function autocorrelation(M, v, a, ns)
    correlations = Float64[];
    # we can take vector n values as input
    for n in ns
        # Sum over M-n terms (v[i+n,a]-v̂[a])(v[i,a]-v̂[a]), divide by (M-n),
        #   is the same as taking the mean of those terms for i in [1:M-n].
        push!(correlations, mean(
            ( v[(n+1):M,a] - mean(v[:,a]) ) .* ( v[1:(M-n),a] - mean(v[:,a]) )
        ));
    end
    return correlations;
end

# find autocorrelation time for datasets v, M samples, cutoff ncut
function τ(M, v, ncut);
    ns = [0:ncut];

    # Find autocorrelation times
    tau = Float64[];
    ndims(v) == 1 ? ndatasets = 1 : ndatasets = size(v)[2];
    for a in [1:ndatasets]
        # Just sum over the normalized autocorrelations
        corra = autocorrelation(M, v, a, ns);
        corra /= corra[1];                 # Normalize
        push!(tau,sum(corra) - 0.5);       # Get autocorrelation time
    end
    return tau;
end

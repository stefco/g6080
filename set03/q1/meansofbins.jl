# Calculate sample means for each bin and dataset
function meansofbins(N, M, v)
    cuts = [N:N:M];         # Cut points for each bucket
    means = Float64[];      # Array holding mean values
    start = 1;
    for cut in cuts
        push!(means,mean(v[start:cut]));
        start = cut + 1;
    end
    return means;
end

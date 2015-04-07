# jacknife sample some set
function jackniferesample(va)
    return ( sum(va) .- va ) / (length(va) - 1);
end

# get variance of set using jacknife resample
function jacknifeσ(va::Array{Float64, 1})
    va′ = jackniferesample(va);
    x = va′ - mean( va′ );
    return sqrt( dot(x,x) * (1 - 1/length(va′)) );
end

# Define the jnestimator struct
type jnestimator
    v̄::Array{Float64,2}
    v̄′::Array{Float64,2}
    μ::Array{Float64,1}
    σ::Array{Float64,1}
end

# Function taking a dataset, sample size, and bin size and
#   returning a struct with the means of the bins (v̄), the
#   jacknife resampled means (v̄′), the mean of the sample (μ),
#   and the estimated standard deviation of the dataset (σ).
function jacknife(v, N, b)
    # Handle at least two datasets at a time (for now)
    nsets = size(v)[2];

    # Split the sample of size N into chunks of size b
    nbins = int(floor(N/b));
    include("../q1/meansofbins.jl");
    
    # Jacknife estimate the averages μ and std dev σ.
    v̄ = zeros(nbins,nsets);
    v̄′ = zeros(nbins,nsets);
    μ = zeros(nsets);
    σ = zeros(nsets);
    for i in [1:nsets]
        v̄[:,i] = meansofbins(b, N, v[:,i]);
        v̄′[:,i] = jackniferesample(v̄[:,i]);        # resample 
        μ[i] = mean(v̄′[:,i]);                      # estimate v̄
        σ[i] = jacknifeσ(v̄[:,i]);                     # estimate σ
    end
    return jnestimator(v̄,v̄′,μ,σ);
end

# estimate σ for vectorized f(v̄) using jacknife estimator input
function jacknifeσ(f::Function, j::jnestimator)
    x = f(j.v̄′) .- f(transpose(j.μ));
    return sqrt( dot(x,x) * ( 1 - 1/length(x) ) );
end


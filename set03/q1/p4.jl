## PART 4
#  Find the integrated autocorrelation times

include("autocorrelation.jl");

println("\nPart 4");
println("\n\tPick ncut = 200 and calculate autocorrelation times");

# Based on result of part 3
ncut = 150;
ns = [0:ncut];

# Find autocorrelation times
tau = Float64[];
for a in [1:5]
  # Just sum over the normalized autocorrelations
  corra = autocorrelation(M, v, a, ns);
  corra /= corra[1];                          # Normalize
  push!(tau,sum(corra) - 0.5);                # Get autocorrelation time
end

println("\n\tAutocorrelation times:\n");
for i in [1:5]
  println("\t\tτ̂v$i:\t$(tau[i])");
end


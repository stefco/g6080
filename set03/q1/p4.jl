## PART 4
#  Find the integrated autocorrelation times

include("autocorrelation.jl");

println("\nPart 4");
println("\n\tPick ncut = 200 and calculate autocorrelation times");

# Based on result of part 3
ncut = 150;

tau = τ(M, v, ncut);

println("\n\tAutocorrelation times:\n");
for i in [1:5]
  println("\t\tτ̂v$i:\t$(tau[i])");
end


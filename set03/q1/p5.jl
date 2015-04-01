## PART 5

println("\nPart 5");
println("\n\tCalculating σ̂ for each a");

# Calculate squared standard deviations
σsq = Float64[];
σ = Float64[];
for a in [1:5]
  σsqa = dot( v[:,a] .- v̄̂[a] , v[:,a] .- v̄̂[a] ) / (M-1);
  σa = sqrt(σsqa);
  println("\t\tσ̂[$a] = $σa");
  push!(σ, σa);
  push!(σsq, σsqa);
end

println("\n\tRecall, from Part 2:");
println("\n\tN1 = $N1\n\tN2 = $N2\n");
println("\n\tWe expect:");
println("\n\t\t2 * τ̂[a] * ( σ̂[a]^2 / σ̂1[a]^2 ) = N1 = $N1")
println("\n\tand")
println("\n\t\t2 * τ̂[a] * ( σ̂[a]^2 / σ̂2[a]^2 ) = N2 = $N2")

println("\n\tTest that hypothesis...");
println("\n\tfor N1:");

for a in [1:5]
  expectn = 2 * tau[a] * σsq[a] / σ̂1[a]^2;
  println("\n\t\t2 * τ̂[a] * ( σ̂[a]^2 / σ̂1[a]^2 ) = N1 = $N1");
end

println("\n\t...and for N2:");
for a in [1:5]
  expectn = 2 * tau[a] * σsq[a] / σ̂2[a]^2;
  println("\n\t\t2 * τ̂[a] * (  σ̂[a]^2 / σ̂2[a]^2 ) = N2 = $N2");
end

println("\n\tWhich is exactly what we would expect.");
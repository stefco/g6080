## PART 6

println("\nPart 6");
# Calculate the true covariance matrix for the data

println("\n\tCalculating the true covariance matrix, ĉ:");
D = v - ones(M) * transpose(v̄̂);
ĉ = transpose(D) * D ./ M;

println("\n\tCalculating the normalized covariance matrix, ρ̂:");
ρ̂ = ĉ ./ ( σ * transpose(σ) )
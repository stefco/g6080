## PART 7

include("autocorrelation.jl");

println("\nPart 7");
println("\n\tPick subjets v1 and v2 of each dataset,");
println("\tsizes N1 = 1,000 & N2 = 10,000.");

v1 = v[[1:N1],:];
v2 = v[[1:N2],:];

println("\n\tEstimate autocorrelations");
println("\n\t...for N1...");
xcorr = [0:10:400];
ycorr1 = zeros( length( xcorr ), 5 );
for a in [1:5]
  # Calculate autocorrelations
  ycorr1[:,a] = autocorrelation(N1, v1, a, xcorr);
  # Normalize by dividing through by Ĉ[1,a]
  ycorr1[:,a] /= ycorr1[1,a];
end
println("\n\t...for N2 (See plots below)");
ycorr2 = zeros( length( xcorr ), 5 );
for a in [1:5]
  # Calculate autocorrelations
  ycorr2[:,a] = autocorrelation(N2, v2, a, xcorr);
  # Normalize by dividing through by Ĉ[1,a]
  ycorr2[:,a] /= ycorr2[1,a];
end

# Make plots for each alpha
println("\n\tPlacing estimated autocorrelation plots 'plots1' and 'plots2':");
using Gadfly
plots1 = Plot[];
plots2 = Plot[];
for i in [1:5]
  # Plot
  push!(plots1, Gadfly.plot ( x = xcorr, y = ycorr1[:,i], Geom.line, Geom.point,
      Guide.xlabel("Separation"), Guide.ylabel("Normalized Autocorrelation"),
      Guide.title(join(["Estimated Autocorrelation for Dataset",string(i),
                        "with Sample Size N1 = $N1"]," "," "))))
  push!(plots2, Gadfly.plot ( x = xcorr, y = ycorr2[:,i], Geom.line, Geom.point,
      Guide.xlabel("Separation"), Guide.ylabel("Normalized Autocorrelation"),
      Guide.title(join(["Estimated Autocorrelation for Dataset",string(i),
                        "with Sample Size N2 = $N2"]," "," "))))
end

# Estimate autocorrelation times
println("\n\tEstimating autocorrelation times τ1 and τ2");
τ1 = Float64[];
τ2 = Float64[];
for a in [1:5]
  corra = autocorrelation(N1, v1, a, ns);
  corra /= corra[1];                          # Normalize
  push!(τ1,sum(corra) - 0.5);                 # Get autocorrelation time
end
for a in [1:5]
  corra = autocorrelation(N2, v2, a, ns);
  corra /= corra[1];                          # Normalize
  push!(τ2,sum(corra) - 0.5);                 # Get autocorrelation time
end

println("\n\tEstimated Autocorrelation Times for N1 = $N1:\n");
for i in [1:5]
  println("\t\tτ̂1v$i:\t$(τ1[i])");
end
println("\n\tEstimated Autocorrelation Times for N2 = $N2:\n");
for i in [1:5]
  println("\t\tτ̂2v$i:\t$(τ2[i])");
end
println("\n\tThe estimated autocorrelation times are clearly far better for the larger sample.");

# Estimate standard deviation
σ1 = Float64[];
σ2 = Float64[];
println("\n\tCreating estimators σ1 for the standard deviation for N1 = $N1");
for a in [1:5]
  diff1 = v1[:,a] .- mean(v1[:,a]);
  σ1a = sqrt( dot( diff1, diff1 ) / (N1-1) );
  println("\t\tσ̂1[$a] = $σ1a");
  push!(σ1, σ1a);
end
println("\n\tCreating estimators σ2 for the standard deviation for N2 = $N2");
for a in [1:5]
  diff2 = v2[:,a] .- mean(v2[:,a]);
  σ2a = sqrt( dot( diff2, diff2 ) / (N2-1) );
  println("\t\tσ̂2[$a] = $σ2a");
  push!(σ2, σ2a);
end

# Estimate standard deviation of the sample mean
σN1 = sqrt( 2 / N1 .* τ1 ) .* σ1;
σN2 = sqrt( 2 / N2 .* τ2 ) .* σ2;
println("\n\tEstimating standard deviations of the sample mean for N1 = $N1:");
for i in [1:5]
  println("\n\t\tσ̂N1[$i] = $(σN1[i])");
end
println("\n\tEstimating standard deviations of the sample mean for N2 = $N2:");
for i in [1:5]
  println("\n\t\tσ̂N2[$i] = $(σN2[i])");
end

# Compare to the actual measured values of the standard deviations of the sample mean
sigmean1 = σ̂1./σN1;
sigmean2 = σ̂2./σN2;
println("\n\tTake ratio of measured standard deviation of sample means over the estimated");
println("\tstandard deviation just found in order to measure the accuracy of the estimate.");
println("\n\tFind the ratio for N1 = $N1:")
for i in [1:5]
  println("\n\t\tσ̂1[$i] / σ̂N1[$i] = $(sigmean1[i])");
end
println("\n\tFind the ratio for N2 = $N2:")
for i in [1:5]
  println("\n\t\tσ̂2[$i] / σ̂N2[$i] = $(sigmean2[i])");
end

println("\n\tFor both sample sizes, these estimates are close to the true values of");
println("\tthe standard deviation of the sample means, though the agreement is far");
println("\tbetter for the larger sample size.");

# Find sample means
v̄1 = zeros(5);
v̄2 = zeros(5);
for i in [1:5]
  v̄1[i] = mean(v1[:,i]);
  v̄2[i] = mean(v2[:,i]);
end

# Finally, compare the normalized covariance matrices ρ1 and ρ2 to the true value ρ̂.
println("\n\tFinally, estimate the normalized covariance matrices ρ1 and ρ2.\n");
println("\tEstimating the covariance matrix c1 for sample size N1 = $N1");
D1 = v1 - ones(N1) * transpose(v̄1);
c1 = transpose(D1) * D1 / N1;
println("\tEstimating the covariance matrix c2 for sample size N2 = $N2");
D2 = v2 - ones(N2) * transpose(v̄2);
c2 = transpose(D2) * D2 / N2;

println("\tEstimating the normalized covariant matrix ρ1 for N1 = $N1");
ρ1 = c1 ./ ( σ1 * transpose(σ1) );
println("\tEstimating the normalized covariant matrix ρ2 for N2 = $N2");
ρ2 = c2 ./ ( σ2 * transpose(σ2) );

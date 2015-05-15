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


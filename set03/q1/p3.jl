## PART 3

include("autocorrelation.jl");

println("\nPart 3");
println("\n\tCalculating autocorrelations");

# Create x and y arrays for plotting
xcorr = [0:10:400];
ycorr = zeros( length( xcorr ), 5 );
for a in [1:5]
  # Calculate autocorrelations
  ycorr[:,a] = autocorrelation(M, v, a, xcorr);
  # Normalize by dividing through by Ĉ[1,a]
  ycorr[:,a] /= ycorr[1,a];
end

# Make plots for each alpha
println("\n\tMaking autocorrelation plots and placing in 'plots' vector");
using Gadfly
plots = Plot[];
for i in [1:5]
  # Plot
  push!(plots, Gadfly.plot ( x = xcorr, y = ycorr[:,i], Geom.line, Geom.point,
      Guide.xlabel("Separation"), Guide.ylabel("Autocorrelation"),
      Guide.title(join(["Autocorrelation for Dataset",string(i)]," "))))
end

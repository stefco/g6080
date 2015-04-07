## PART 2

# Find integrated autocorrelation times

include("../q1/autocorrelation.jl");

println("\nPart 2");
println("\n\tFirst, plot autocorrelation to estimate τint");

# Same (almost) as q1p3:
Mmd = length(ts);

# Create x and y arrays for plotting
xcorr = [0:10:1000];
ycorr = zeros( length( xcorr ), 3 );
for a in [1:3]
  # Calculate autocorrelations
  ycorr[:,a] = autocorrelation(Mmd, md, a, xcorr);
  # Normalize by dividing through by Ĉ[1,a]
  ycorr[:,a] /= ycorr[1,a];
end

# Make plots for each alpha
println("\n\tMaking autocorrelation plots and placing in 'mdautocorr' vector");
using Gadfly
mdautocorr = Plot[];
quantity = ["Temperature","Potential Energy","Virial"];
for i in [1:3]
  # Plot
  push!(mdautocorr, Gadfly.plot ( x = xcorr, y = ycorr[:,i], Geom.line, Geom.point,
      Guide.xlabel("Step Separation"), Guide.ylabel("Autocorrelation"),
      Guide.title(join(["Autocorrelation for ",quantity[i]]," "))))
end

# Find integrated autocorrelation times τint (similar to q1p4)

# Pick ncut for temperature, potential energy, and the virial
#   based on observed graph of autocorrelation
ncutt = 150;
ncutp = 150;
ncutv = 150;

τtemperature = τ(Mmd, ts, ncutt);
τpressure = τ(Mmd, ps, ncutp);
τvirial = τ(Mmd, vs, ncutv);

τmd = [τtemperature, τpressure, τvirial];

println("\n\tAutocorrelation times:\n");
for i in [1:3]
  println("\t\tτ̂v",quantity[i]":\t",τmd[i]);
end


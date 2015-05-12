## PART 2 -- Note: you have to look at p1.jl for this to make sense

# Find integrated autocorrelation times

include("../q1/autocorrelation.jl");

println("\nPart 2");
println("\n\tFirst, plot autocorrelation to estimate τint");

# Fuse md trials
md = [md_1069 md_1304]

# Number of quantities to look at
N = 6

# Same (almost) as q1p3:
lengthmd = size(md)[1];

# Create x and y arrays for plotting
xcorr = [0:10:1000];
ycorr = zeros( length( xcorr ), N );
for a in [1:N]
  # Calculate autocorrelations
  ycorr[:,a] = autocorrelation(lengthmd, md, a, xcorr);
  # Normalize by dividing through by Ĉ[1,a]
  ycorr[:,a] /= ycorr[1,a];
end

# Make plots for each alpha
println("\n\tMaking autocorrelation plots and placing in 'mdautocorr' vector");
using Gadfly
mdautocorr = Plot[];
quantity = ["Temperature at T = 1.069", "Potential Energy at T = 1.069",
            "Virial at T = 1.069", "Temperature at T = 1.304",
            "Potential Energy at T = 1.304", "Virial at T = 1.304"];

for i in [1:N]
  # Plot
  push!(mdautocorr, Gadfly.plot ( x = xcorr, y = ycorr[:,i], Geom.line, Geom.point,
      Guide.xlabel("Step Separation"), Guide.ylabel("Autocorrelation"),
      Guide.title(join(["Autocorrelation for ",quantity[i]]," "))))
end

# Find integrated autocorrelation times τint (similar to q1p4)

# Pick ncut for temperature, potential energy, and the virial
#   based on observed graph of autocorrelation
ncutt_1069 = 150;
ncutp_1069 = 150;
ncutv_1069 = 150;
ncutt_1304 = 150;
ncutp_1304 = 150;
ncutv_1304 = 150;

τtemp_1069 = τ(lengthmd, temp_1069, ncutt_1069);
τpe_1069 = τ(lengthmd, pe_1069, ncutp_1069);
τvir_1069 = τ(lengthmd, vir_1069, ncutv_1069);
τtemp_1304 = τ(lengthmd, temp_1304, ncutt_1304);
τpe_1304 = τ(lengthmd, pe_1304, ncutp_1304);
τvir_1304 = τ(lengthmd, vir_1304, ncutv_1304);

τmd = [τtemp_1069, τpe_1069, τvir_1069, τtemp_1304, τpe_1304, τvir_1304]

println("\n\tAutocorrelation times:\n");
for i in [1:N]
  println("\t\tτ̂v",quantity[i]":\t",τmd[i]);
end


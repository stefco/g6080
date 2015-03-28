## PART 2 

include("meansofbins.jl");

# Sample sizes
N1 = 1000;       # One thousand
N2 = 10000;      # Ten thousand
println("\nPart 2");
println("\n\tN1 = $N1\n\tN2 = $N2\n")

# Leave space for the means v̄1 and v̄2
v̄1 = zeros(int64(M/N1),5);
v̄2 = zeros(int64(M/N2),5);

# Save all those means
for a in [1:5]
    v̄1[:,a] = meansofbins(N1, M, v[:,a]);
    v̄2[:,a] = meansofbins(N2, M, v[:,a]);
end

# Declare standard devs of the sample means
σ̂1 = Float64[];
σ̂2 = Float64[];

# Calculate standard devs for sample means
for a in [1:5]
    push!(σ̂1,std(v̄1[:,a]));
    push!(σ̂2,std(v̄2[:,a]));
end

# Print the results
println("\t Standard deviations of the sample means are:\n");
println("\tσ̂1 = $σ̂1");
println("\tσ̂2 = $σ̂2");

# See if this is what we'd expect
println("\n\tWe expect (σ̂1./σ̂2).^2 = N2/N1 = 10:\n");
for i in [1:5]
    ratiosqr = (σ̂1[i]/σ̂2[i])^2
    println("\t(σ̂1[$i]/σ̂2[$i])^2 = $ratiosqr");
end
println("\n\t...which is pretty close.");

println("\nSee below for histograms.");

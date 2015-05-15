## PART 4

println("\nPart 4");
include("../q2/jacknife.jl");

# Set bin size b to the largest autocorrelation time for each temp
numBodies = 1000;
println("\n\tFinding number of steps between independent data recordings\n")
@show bp_1069 = 2maximum([τtemp_1069, τpe_1069, τvir_1069])
@show bp_1304 = 2maximum([τtemp_1304, τpe_1304, τvir_1304])

b=10
println("\n\tPick a bin size of $b\n")

# Function for calculating pressure at various temperatures
function pressure(numBodies, temperatures, virials)
    # Take the time average of the virial
    virialavg = mean(virials);

    # Set Boltzmann constant
    kB = 1.3806e-23;

    kBTs = temperatures * kB;

    return 1 - virialavg * ( 1 ./ (6 * numBodies * kBTs ) );
end

println("\n\tEstimating σpressure:");

mdjack_1069 = jacknife(md_1069, lengthmd, b);
mdjack_1304 = jacknife(md_1304, lengthmd, b);

# Define a function that we can pass to jacknifeσ
pressureN(md) = pressure(numBodies, md[:,1], md[:,3]);

# Estimate σpressure
σpressure_1069 = jacknifeσ(pressureN, mdjack_1069);
σpressure_1304 = jacknifeσ(pressureN, mdjack_1304);

println("\n\tAverage pressure for T = 1.069 is ", mean(pressureN(md_1069)))
println("\n\tAverage pressure for T = 1.304 is ", mean(pressureN(md_1304)))

println("\n\tStandard deviation of pressure for T = 1.069: $σpressure_1069")
println("\n\tStandard deviation of pressure for T = 1.304: $σpressure_1304")

# Should also plot pressure

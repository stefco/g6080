## PART 4

println("\nPart 4");
include("../q2/jacknife.jl");

# Set bin size b to the largest autocorrelation time for each temp
numBodies = 1000;
bp_1069 = 2maximum(ncutt_1069, ncutp_1069, ncutv_1069)
bp_1304 = 2maximum(ncutt_1304, ncutp_1304, ncutv_1304)

# Function for calculating pressure at various temperatures
function pressure(numBodies, temperatures, virials)
    # Take the time average of the virial
    virialavg = mean(vs);

    # Set Boltzmann constant
    kB = 1.3806e-23;

    kBTs = temperatures * kB;

    return 1 - virialavg * ( 1 ./ (6 * numBodies * kBTs ) );
end

println("\n\tEstimating σpressure for bin size b=$b:");

mdjack_1069 = jacknife(md_1069, lengthmd, bp_1069);
mdjack_1304 = jacknife(md_1304, lengthmd, bp_1304);

# Define a function that we can pass to jacknifeσ
pressureN(md) = pressure(numBodies, md[:,1], md[:,3]);

# Estimate σpressure
σpressure_1069 = jacknifeσ(pressureN, mdjack[1:3]);
σpressure_1304 = jacknifeσ(pressureN, mdjack[4:6]);

println("\n\tStandard deviation of pressure for T = 1.069: $σpressure_1069")
println("\n\tStandard deviation of pressure for T = 1.304: $σpressure_1304")

# Should also plot pressure

## PART 4

println("\nPart 4");
include("../q2/jacknife.jl");

# Set bin size b to the largest autocorrelation time
numBodies = 1000;
bp = 150;

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

mdjack = jacknife(md, Mmd, bp);

# Define a function that we can pass to jacknifeσ
pressureN(md) = pressure(numBodies, md[:,1], md[:,3]);

# Estimate σpressure
σpressure = jacknifeσ(pressureN, mdjack);

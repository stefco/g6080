## PART 1

println("~~ QUESTION 4 ~~");
println("\nPart 1");

# Read in values from previous simulation (but only after thermalization):
println("\n\tReading in values from MD simulation...");

temps_1069 = readdlm("q4/temps-rho-0.75-temp-1.069.dat")[401:end]
temps_1304 = readdlm("q4/temps-rho-0.75-temp-1.304.dat")[401:end]

pe_1069 = readdlm("q4/pe-rho-0.75-temp-1.069.dat")[401:end]
pe_1304 = readdlm("q4/pe-rho-0.75-temp-1.304.dat")[401:end]

vir_1069 = readdlm("q4/vir-rho-0.75-temp-1.069.dat")[401:end]
vir_1304 = readdlm("q4/vir-rho-0.75-temp-1.304.dat")[401:end]

# Slap it all into a single array so that we can reuuse
#   old code
md_1069 = [temps_1069 pe_1069 vir_1069];
md_1304 = [temps_1304 pe_1304 vir_1304];

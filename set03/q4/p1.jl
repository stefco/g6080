## PART 1

println("~~ QUESTION 4 ~~");
println("\nPart 1");

# Read in values from previous simulation:
println("\n\tReading in values from MD simulation...");
ts = readdlm("q4/ts");
ps = readdlm("q4/ps");
vs = readdlm("q4/vs");

# Slap it all into a single array so that we can reuuse
#   old code
md = [ts ps vs];

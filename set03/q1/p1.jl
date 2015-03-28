## PART 1

# Declare function for loading sample data
include("readarray.jl");

# Set the length of each array
M = 1600000;

# Read in values of each array
v1 = readarray("v1");
v2 = readarray("v2");
v3 = readarray("v3");
v4 = readarray("v4");
v5 = readarray("v5");
v = [v1 v2 v3 v4 v5];

# Take mean values of each array
v̄̂ = Float64[];
push!(v̄̂,mean(v[:,1]));
push!(v̄̂,mean(v[:,2]));
push!(v̄̂,mean(v[:,3]));
push!(v̄̂,mean(v[:,4]));
push!(v̄̂,mean(v[:,5]));

# Print the result nicely
println("\n~~ QUESTION 1 ~~\n");
println("Part 1");
println("\tThe true means v̄̂1, v̄̂2, ..., v̄̂5 are:\n");
println("\tv̄̂[1] = $(v̄̂[1])");
println("\tv̄̂[2] = $(v̄̂[2])");
println("\tv̄̂[3] = $(v̄̂[3])");
println("\tv̄̂[4] = $(v̄̂[4])");
println("\tv̄̂[5] = $(v̄̂[5])");


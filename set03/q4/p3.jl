## PART 3

println("\nPart 3");

# Estimate the covariance matrix for these quantities
println("\n\tCalculating the covariance matrices, Cmoldy_1069 and Cmoldy_1304");

# Calculate the mean values of each quantity
mdavg_1069 = zeros(3);
mdavg_1304 = zeros(3);
for n in 1:3
    mdavg_1069[n] = mean(md_1069[:,n]);
    mdavg_1304[n] = mean(md_1304[:,n]);
end

# Calculate the covariance for T = 1.069...
println("\n\tCovariance matrix for T = 1.069:\n")
D = md_1069 - ones(lengthmd) * transpose(mdavg_1069);
Cmoldy_1069 = transpose(D) * D ./ lengthmd
show(Cmoldy_1069)

# ...and for T = 1.304.
println("\n\n\tCovariance matrix for T = 1.304\n")
D = md_1304 - ones(lengthmd) * transpose(mdavg_1304);
Cmoldy_1304 = transpose(D) * D ./ lengthmd;
show(Cmoldy_1304)

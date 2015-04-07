## PART 3

println("\nPart 3");

# Estimate the covariance matrix for these quantities
println("\n\tCalculating the covariance matrix, Cmoldy:");

# Calculate the mean values of each quantity
mdavg = zeros(3);
mdavg[1] = mean(md[:,1]);
mdavg[2] = mean(md[:,2]);
mdavg[3] = mean(md[:,3]);

# Calculate the covariance
D = md - ones(Mmd) * transpose(mdavg);
Cmoldy = transpose(D) * D ./ Mmd;

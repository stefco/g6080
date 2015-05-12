## PART 3

println("\nPart 3");

# Estimate the covariance matrix for these quantities
println("\n\tCalculating the covariance matrix, Cmoldy:");

# Calculate the mean values of each quantity
mdavg = zeros(N);
for n in 1:N
    mdavg[n] = mean(md[:,n]);
end

# Calculate the covariance
D = md - ones(lengthmd) * transpose(mdavg);
Cmoldy = transpose(D) * D ./ lengthmd;

## PART 1

include("partition.jl");
include("../q1/autocorrelation.jl");

println("\n\tPick ncut = 150, as before");
ncut = 150;

# Partition the set into M/N bins
println("\n\tPartitioning v into M/N bins for N1 = $N1 and N2 = $N2");
vpartsN1 = partition(v, M, N1);
vpartsN2 = partition(v, M, N2);

println("\n\tCounting number of bins");
npartsN1 = size(vpartsN1)[3];
npartsN2 = size(vpartsN2)[3];

println("\n\tCalculating τint for M/N1 bins");
τintN1 = zeros(npartsN1,5);
for npart in [1:npartsN1]
    vpart = vpartsN1[:,:,npart];
    τintN1[npart,:] = τ(N1, vpart, ncut);
end

println("\n\tCalculating στint for M/N1 bins");
μτintN1 = zeros(5);
στintN1 = zeros(5);
for i in [1:5]
    μτintN1[i] = mean(τintN1[:,i]);
    στintN1[i] = std(τintN1[:,i]);
end

println("\n\tCalculating τint for M/N2 bins");
τintN2 = zeros(npartsN2,5);
for npart in [1:npartsN2]
    vpart = vpartsN2[:,:,npart];
    τintN2[npart,:] = τ(N2, vpart, ncut);
end

println("\n\tCalculating στint for M/N2 bins");
μτintN2 = zeros(5);
στintN2 = zeros(5);
for i in [1:5]
    μτintN2[i] = mean(τintN2[:,i]);
    στintN2[i] = std(τintN2[:,i]);
end

## PART 3

println("Part 3");
println("\n\tEstimating mean and standard deviation using"); 
println("\tjacknife resampling with bin size of 40.");

include("jacknife.jl");

b = 40;
jest = jacknife(v, N2, b);

println("\n\tEstimated means, using N = $N2 and b = $b:\n");
for i in [1:5]
    println("\t\tdataset $i: $(jest.μ[i])");
end

println("\n\tEstimated standard deviations:\n");
for i in [1:5]
    println("\t\tdataset $i: $(jest.σ[i])");
end

println("\n\tCalculating dependence of jacknife σ estimator on b.");
bvals = [1:100];
σb = zeros(length(bvals),5);
for i in [1:length(bvals)]
    jestb = jacknife(v, N2, bvals[i]);
    σb[i,:] = jestb.σ;
end


println("\n\tGenerating plots in σvsbplots to show that dependence.");
using Gadfly
σvsbplots = Plot[];
ttl = "Estimated Variance vs. Bin Size for Dataset";
for i in [1:5]
    push!(σvsbplots, Gadfly.plot(
        x = bvals,
        y = σb[:,i],
        Geom.line,
        Geom.point,
        Guide.xlabel("Bin size b"),
        Guide.ylabel("Estimated Variance"),
        Guide.title(join([ttl,string(i)]," "))
    ));
end


## PART 4

println("\nPart 4");

bvals = [b-3:b+3];
println("\n\tEstimating σ for f1..f3 with bin around b=$b");

σf1jack = Float64[];
σf2jack = Float64[];
σf3jack = Float64[];

for bb in bvals
    jestb = jacknife(v, N2, bb);
    push!(σf1jack, jacknifeσ(f1, jestb));
    push!(σf2jack, jacknifeσ(f2, jestb));
    push!(σf3jack, jacknifeσ(f3, jestb));
end

σfjack = cat(2, σf1jack, σf2jack, σf3jack);

println("\n\tGenerating plots in σfvsb to show σf dependence on b.");
using Gadfly
σfvsb = Plot[];
for i in [1:3]
    push!(σfvsb, Gadfly.plot(
        x = bvals,
        y = σfjack[:,i],
        Geom.line,
        Geom.point,
        Guide.xlabel("Bin size b"),
        Guide.ylabel("Estimated Variance"),
        Guide.title("Estimated Variance of function f$i vs. Bin Size")
    ));
end

println("\n\tCompare σfjack to naive σ̂ from part 1:\n")
println("\t\tσ̂f1N2/σfjack[4,1] = $(σ̂f1N2/σfjack[4,1])");
println("\t\tσ̂f2N2/σfjack[4,2] = $(σ̂f2N2/σfjack[4,2])");
println("\t\tσ̂f3N2/σfjack[4,3] = $(σ̂f3N2/σfjack[4,3])");

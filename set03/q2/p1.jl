## PART 1

println("Part 1");

# Define functions for this part
function f1(v̄)
    return v̄[:,1] ./ v̄[:,2];
end

function f2(v̄)
    return exp( v̄[:,3] - v̄[:,4] );
end

function f3(v̄)
    logs = Float64[];
    for vv in v̄[:,5]
        push!(logs,log(vv));
    end
    return ( v̄[:,1]./v̄[:,2] + v̄[:,3]./v̄[:,4] ) .* logs;
end

# Calculate f1, f2, and f3 for these means
println("RESUMING QUESTION 2");
println("\n\tCalculating functions of means f1v1, f2v1, ...");
f1N1 = f1(v̄1);
f2N1 = f2(v̄1);
f3N1 = f3(v̄1);

f1N2 = f1(v̄2);
f2N2 = f2(v̄2);
f3N2 = f3(v̄2);

# Find standard deviations
println("\n\tCalculating standard deviations of the functions");
σ̂truef1N1 = std(f1N1);
σ̂truef2N1 = std(f2N1);
σ̂truef3N1 = std(f3N1);

σ̂truef1N2 = std(f1N2);
σ̂truef2N2 = std(f2N2);
σ̂truef3N2 = std(f3N2);


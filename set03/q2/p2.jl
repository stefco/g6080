## PART 2

# Define the error functions

function σf1(v̄̂,σ̂)
    ddvs = zeros(5);
    ddvs[1] = 1 / v̄̂[2]^2;
    ddvs[2] = v̄̂[1]^2 / v̄̂[2]^4;
    return sqrt( dot( ddvs, σ̂.^2) );
end

function σf2(v̄̂,σ̂)
    ddvs = zeros(5);
    ddvs[3] = exp(v̄̂[3] - v̄̂[4])^2;
    ddvs[4] = exp(v̄̂[3] - v̄̂[4])^2;
    return sqrt( dot( ddvs, σ̂.^2) );
end

function σf3(v̄̂,σ̂)
    lv5s = log(v̄̂[5])^2;
    ddvs = zeros(5);
    ddvs[1] = lv5s / v̄̂[2]^2;
    ddvs[2] = lv5s * v̄̂[1]^2 / v̄̂[2]^4;
    ddvs[3] = lv5s / v̄̂[4]^2;
    ddvs[4] = lv5s * v̄̂[3]^2 / v̄̂[4]^4;
    ddvs[5] = ( ( v̄̂[1]/v̄̂[2] + v̄̂[3]/v̄̂[4] ) / v̄̂[5] )^2;
    return sqrt( dot( ddvs, σ̂.^2) );
end

println("Calculating σ̂f1, σ̂f2, σ̂f3 for sample sizes N1, N2");
σ̂f1N1 = σf1(v̄̂,σ̂1); 
σ̂f2N1 = σf2(v̄̂,σ̂1);
σ̂f3N1 = σf3(v̄̂,σ̂1);

σ̂f1N2 = σf1(v̄̂,σ̂2); 
σ̂f2N2 = σf2(v̄̂,σ̂2);
σ̂f3N2 = σf3(v̄̂,σ̂2);

println("\n\tNaive standard deviations for sample size 1,000");
println("\tσ̂f1N1 = $σ̂f1N1"); 
println("\tσ̂f2N1 = $σ̂f2N1");
println("\tσ̂f3N1 = $σ̂f3N1");

println("\n\tNaive standard deviations for sample size 10,000");
println("\tσ̂f1N2 = $σ̂f1N2"); 
println("\tσ̂f2N2 = $σ̂f2N2");
println("\tσ̂f3N2 = $σ̂f3N2");

println("\n\tTrue standard deviations for sample size 1,000");
println("\tσ̂truef1N1 = $σ̂truef1N1"); 
println("\tσ̂truef2N1 = $σ̂truef2N1");
println("\tσ̂truef3N1 = $σ̂truef3N1");

println("\n\tTrue standard deviations for sample size 10,000");
println("\tσ̂truef1N2 = $σ̂truef1N2"); 
println("\tσ̂truef2N2 = $σ̂truef2N2");
println("\tσ̂truef3N2 = $σ̂truef3N2");


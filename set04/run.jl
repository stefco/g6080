include("configure.jl")

σs = Array(Float64,length(Ts),length(Ns))   # standard devs of magnetizations
absμs = Array(Float64,length(Ts),length(Ns))   # means of abs val of magnetizations
μs = Array(Float64,length(Ts),length(Ns))   # means of abs val of magnetizations
acc = Array(Float64,length(Ts),length(Ns))  # accept ratio for these parameters

mags = Array(Float64, steps)
absmags = Array(Float64, steps)

for n in 1:length(Ns)
    for b in 2:(length(Bs) - 1)
        for t in 1:length(Ts)
            @time mags, acc[t,n] = cluster(Ns[n], Ns[n], steps, J, Bs[b], Ts[t]);
    
            absmags = abs(mags);
            absμs[t,n] = mean(absmags);
            μs[t,n] = mean(mags);
            σs[t,n] = std(mags);
            
            writedlm("absmeans-$(Bs[b])", absμs)
            writedlm("means-$(Bs[b])", μs)
            writedlm("stddevs-$(Bs[b])", σs)
            writedlm("acceptratios-$(Bs[b])", acc)
        end
    end
end


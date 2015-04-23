include("configure.jl")

σs = Array(Float64,length(Ts),length(Bs))   # standard devs of magnetizations
μs = Array(Float64,length(Ts),length(Bs))   # means of abs val of magnetizations
acc = Array(Float64,length(Ts),length(Bs))  # accept ratio for these parameters

mags = Array(Float64, steps)
absmags = Array(Float64, steps)

for b in 1:length(Bs)
    for t in 1:length(Ts)
        mags, acc[t,b] = cluster(Nx, Ny, steps, J, Bs[b], Ts[t])
        absmags = abs(mags)
        μs[t,b] = mean(absmags)
        σs[t,b] = std(mags)
    end
end

writedlm("means", μs)
writedlm("stddevs", σs)
writedlm("acceptratios", acc)

include("configure.jl")

corrTs = [2.1, 2.2, 2.25, 2.26, 2.265, 2.269]
μs = Array(Float64,length(mcd),length(corrTs))   # means of abs val of magnetizations
σs = Array(Float64,length(mcd),length(corrTs))   # standard devs of magnetizations
estimators = Array(Float64, length(corrTs))
corrSteps = 20000

# find spatial correlators for multiple T values
for t in 1:length(corrTs)

    @time mags, acc, naivecorrelators, correlators =
            cluster(Nx, Ny, corrSteps, J, 0.0, corrTs[t], mcd)

    # need to find μs and σs for the spatial correlators at diff separations
    μs[:,corrTs[t]] = mean(naivecorrelators, 2)
    σs[:,corrTs[t]] = std(naivecorrelators, 2)
    estimators[t] = mean(correlators)

    writedlm("cormeans", μs)
    writedlm("corstddevs", σs)
    writedlm("estimatormeans", estimatormeans)
end

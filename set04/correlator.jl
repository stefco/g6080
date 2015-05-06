include("configure.jl")

corrTs = [2.1, 2.2, 2.25, 2.26, 2.265, 2.269]
μs = Array(Float64,mcd,length(corrTs))   # means of correlator
μas = Array(Float64,mcd,length(corrTs))  # sample means of expected radius
σs = Array(Float64,mcd,length(corrTs))   # standard devs of correlator
σas = Array(Float64,mcd,length(corrTs))  # sample standard devs of expected radius
estimators = Array(Float64, length(corrTs))
corrSteps = 2000

# find spatial correlators for multiple T values
for t in 1:length(corrTs)

    @time mags, acc, naivecorrelators =    # throw in correlators later on
            cluster(Nx, Ny, corrSteps, J, 0.0, corrTs[t], mcd)

    # need to find μs and σs for the spatial correlators at diff separations
    μs[:,t] = mean(naivecorrelators, 2)
    μas[:,t] = mean(aoft(naivecorrelators), 2)
    σs[:,t] = std(naivecorrelators, 2)
    σas[:,t] = std(aoft(naivecorrelators), 2)
    ### estimators[t] = mean(correlators)

    writedlm("cormeans", μs)
    writedlm("corstddevs", σs)
    ### writedlm("estimatormeans", estimatormeans) # DOESNT WORK ANYWAY
end

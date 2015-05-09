## QUESTION 2

include("md.jl");
qstart = loadmdtrial("../../set04/qstart.mdv")

# Pick interaction distance
δ = 0.1
threads = 2
N = 60
M = 10
P = int(N/M)

println("Running ",threads," threads with ",N," steps each and ",P," measurements")

PE = Array(Float64, P, threads)
KE = Array(Float64, P, threads)
V  = Array(Float64, P, threads)

# Place the particles using result from qtherm
qs = Array(MolecularDynamicsTrial,threads)
for i in 1:threads
    qs[i] = addsteps(qstart, N)     # N new steps
end

# Since the metropolis algorithm is stochastic and aphysical, we can start from 
# the thermalized state and just run several simulations in parallel. Set it up
# for an 8 core google compute engine server. Total of 4.8e4 steps.
@sync @parallel for i in 1:threads
    
    q = qs[i]

    # Store potential and kinetic energy as well as the virial every 10 steps
    pe = Float64[]
    ke = Float64[]
    v  = Float64[]

    # Run for 100 steps at a time, saving values along the way
    while !isfinished(q)
        metropolis!(q, M, δ)    # M steps between saves
        push!(pe, q.pet[q.currentStep])
        push!(ke, q.ket[q.currentStep])
        push!(v,  q.vir[q.currentStep])
    end

    # append results to our collection
    PE[:,i] = pe
    KE[:,i] = ke
    V[:,i]  = v
end

# save our results
writedlm("$(pwd())/set04/pe.dat", PE)
writedlm("$(pwd())/set04/ke.dat", KE)
writedlm("$(pwd())/set04/v.dat",  V)

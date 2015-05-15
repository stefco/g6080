## QUESTION 2

include("md.jl");
qstart = loadmdtrial("../../swp/1069_met_start.mdv")

δpe = 0.1                     # interaction distance
δke = 0.1                     # interaction distance
N = 20000                   # simulation steps
M = 10                      # steps between saved measurements
P = int(N/M)

println("Running ",N," steps for ",P," measurements")

q = addsteps(qstart, N)     # N new steps
pe = Float64[]
ke = Float64[]
v  = Float64[]

# Run for 10 steps at a time, saving values along the way
while !isfinished(q)
    metropolis!(q, M, δpe, δke)    # M steps between saves
    push!(pe, q.pet[q.currentStep])
    push!(ke, q.ket[q.currentStep])
    push!(v,  q.vir[q.currentStep])
end

# save our results
writedlm("../../swp/set04/q1/pe.dat", pe)
writedlm("../../swp/set04/q1/ke.dat", ke)
writedlm("../../swp/set04/q1/v.dat",  v)

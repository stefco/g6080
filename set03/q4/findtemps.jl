# given particle quantity and an array of kinetic energies, find the temperature 
# at each step
function temperatures(ke::Array{Float64}, N::Integer)
    return 2ke/3N
end

function temperatures(fname::String, N::Integer)
    ke = readdlmn(fname)
    return temperatures(ke, N)
end

function temperatures(rname::String, sname::String, N::Integer)
    t = temperatures(rname, N)
    writedlm(sname, t)
    return t
end

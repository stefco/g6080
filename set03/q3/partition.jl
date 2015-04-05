# Partition the first M rows of v into M/N sets of length N
function partition(v::Array{Float64,2}, M::Int64, N::Int64)
    npartitions = int(floor(M/N));
    vparts = zeros(N,5,npartitions);   # preallocate results
    steps = [1:npartitions];
    for i in steps
        binstart = N*(i-1) + 1;
        binend = N*i;
        vparts[:,:,i] = v[binstart:binend,:];
    end
    return vparts;
end

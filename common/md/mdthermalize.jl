function thermalize(filepath::String, r::MolecularDynamicsTrial, N::Int64, 
        n::Int64 )
    for i in 1:n
        println("Thermalizing, step $i of $n, with $N verlet steps per relaxation")
        verlet!(r,N)
        save(filepath,r)
        relax!(r)
        save(filepath,r)
    end
end

function thermalize(filepath::String, r::MolecularDynamicsTrial, N::Int64, 
        n::Int64, index::MolecularDynamicsIndex )
    for i in 1:n
        println("Thermalizing, step $i of $n, with $N verlet steps per relaxation")
        verlet!(r,N,index)
        save(join([filepath,r.currentStep]),r)
        relax!(r)
        save(join([filepath,r.currentStep,"r"]),r)
    end
end

module q2

import Evolve

steps = 16000               # number of simulation steps
samplesteps = 50            # steps between measurements
L = 2                       # width of the box
xmin = -1                   # leftmost point in the box
N = 2000                    # number of subintervals
ϵ = L/N                     # width of each subinterval
ks = [100, 200]              # wavenumbers
σs = [0.05, 0.02]            # wavepacket width
x0 = 0.0                    # wavepacket starting position

# Run the trials
function run()
for k in ks
    for σ in σs
        print("Initializing for k = $k, σ = $σ...")
        V(xs) = map(x->0.0, xs)     # zero potential energy
        AV = Evolve.A(V)            # evolution matrix constructor for square well
        ϕV = Evolve.Φ(k, σ, x0)     # initial state vector constructor for Gaussian packet

        A = AV(L, N, xmin)
        ϕ0 = ϕV(L, N, xmin)
        x = (L/N .* [0:N] + xmin)   # x positions

        print("evolving...")
        ϕ = Evolve.evolve(A, ϕ0, steps)
        phis = [ϕ0 ϕ[:,samplesteps:samplesteps:end]]

        println("done. Writing to file.")
        writedlm("phis_im_$k\_$σ.dat", imag(phis))
        writedlm("phis_re_$k\_$σ.dat", real(phis))
    end
end
end

# Run the trials
for k in ks
    for σ in σs
        print("Reverse run for k = $k, σ = $σ...")
        V(xs) = map(x->0.0, xs)     # zero potential energy
        AV = Evolve.A(V)            # evolution matrix constructor for square well
        ϕV = Evolve.Φ(k, σ, x0)     # initial state vector constructor for Gaussian packet

        A = AV(L, N, xmin)
        ϕ0 = ϕV(L, N, xmin)
        x = (L/N .* [0:N] + xmin)   # x positions

        print("evolving...")
        ϕ = Evolve.evolve(A, ϕ0, 2000)
        phis = [ϕ0 ϕ[:,samplesteps:samplesteps:end]]
        ϕr = conj(ϕ[:,end])         # reverse time evolution by taking conjugate
        ϕ = Evolve.evolve(A, ϕr, 2000)
        phis = [phis ϕ[:,samplesteps:samplesteps:end]]

        println("done. Writing to file.")
        writedlm("phis_rev_im_$k\_$σ.dat", imag(phis))
        writedlm("phis_rev_re_$k\_$σ.dat", real(phis))
    end
end

end

module q5

import Evolve

steps = 8000                # number of simulation steps
samplesteps = 50            # steps between measurements
L = 2                       # width of the box
xmin = -1                   # leftmost point in the box
N = 2000                    # number of subintervals
ϵ = L/N                     # width of each subinterval
h = 1e5                     # energy is 12336.8 for w = 0.02
b = 1
s = 0.005                    # wall thickness
w = 0.02         # pocket widths to try
ks = π/w .* [2.4:0.2:3.2]                  # wavenumbers
σ = 0.05                    # wavepacket width
V0s = [0]         # step potentials
x0 = -0.5                   # wavepacket starting position
Ts = Float64[]
Rs = Float64[]

# Run the trials
for k in ks
    for V0 in V0s
        print("Initializing for k = $k... ")
        V(xs) = map(x -> Evolve.pocket(L, h, b, s, w, xmin, x), xs)      # pocket potential
        AV = Evolve.A(V)            # evolution matrix constructor for square well
        ϕV = Evolve.Φ(k, σ, x0)     # initial state vector constructor for Gaussian packet

        A = AV(L, N, xmin)
        ϕ0 = ϕV(L, N, xmin)
        x = (L/N .* [0:N] + xmin)   # x positions

        # Evolve the system and record positions every so many steps
        print("evolving... ")
        test = x -> (abs(x[2]) > 1e-2 || abs(x[end-1]) > 1e-2)
        ϕ = Evolve.evolve(A, ϕ0, steps, test)
        phis = [ϕ0 ϕ[:,samplesteps:samplesteps:end]]

        # calculate transmission and reflection probabilities
        T = mean(abs(ϕ[N/2+1:end,end]).^2)
        R = mean(abs(ϕ[1:N/2,end]).^2)
        I = T + R
        T /= I
        R /= I
        push!(Ts, T)
        push!(Rs, R)
        print("T = $T, R = $R...")
        ### writedlm("trans_step_$V0", [T,R,I])

        println("done after $(size(ϕ)[2]) steps. Writing to file.")
        writedlm("phis_scatter_im_$k.dat", imag(phis))
        writedlm("phis_scatter_re_$k.dat", real(phis))
    end
end

writedlm("scatter1.dat", [ks Ts Rs])
end

module q4

import Evolve

steps = 32000                # number of simulation steps
samplesteps = 50            # steps between measurements
L = 2                       # width of the box
xmin = -1                   # leftmost point in the box
N = 2000                    # number of subintervals
ϵ = L/N                     # width of each subinterval
k = 200                     # wavenumbers
σ = 0.05                    # wavepacket width
h = 1000                      # energy is 123.368 for w = 0.2
b = 1
ss = [0.001, 0.01, 0.1]                    # wall thickness
w = 0.2                     # pocket width
x0 = -0.5                   # wavepacket starting position

# Functions for square well states and 1st excited state
Φground = x -> √(2/w)*cos(π*x/w)
Φexcited = x -> √(2/w)*sin(2*π*x/w)
Φs = [:Φground, :Φexcited]

# Run the trials
@sync @parallel for Φ0 in Φs
    for s in ss
        print("Initializing for Φ0 = $(string(Φ0)), s = $s... ")
        V(xs) = map(x -> Evolve.pocket(L, h, b, s, w, xmin, x), xs)      # pocket potential
        AV = Evolve.A(V)            # evolution matrix constructor for square well

        A = AV(L, N, xmin)
        Φ = eval(Φ0)
        ϕ0 = zeros(Complex, N+1)
        xpocket = [(0.5L-0.5w):ϵ:(0.5L+0.5w)]
        ϕ0[int((0.5L-0.5w)N/L):int((0.5L+0.5w)N/L)] = Φ(xpocket)
        x = (L/N .* [0:N] + xmin)   # x positions

        # Evolve the system and record positions every so many steps
        print("evolving... ")
        test = x -> (abs(x[2]) > 1e-2 || abs(x[end-1]) > 1e-2)
        ϕ = Evolve.evolve(A, ϕ0, steps, test)
        phis = [ϕ0 ϕ[:,samplesteps:samplesteps:end]]

        # calculate transmission and reflection probabilities
        print("calculating transmission probabilities...")
        Ts = Float64[]
        Rs = Float64[]
        for i in 1:size(phis)[2]
            phi = phis[:,i]
            T = mean(abs(phi[1:int((0.5L-0.5w-s)N/L),end]).^2)
            T += mean(abs(phi[int((0.5L+0.5w+s)N/L):end,end]).^2)
            R = mean(abs(phi[int((0.5L-0.5w)N/L):int((0.5L+0.5w)N/L),end]).^2)
            I = T + R
            T /= I
            R /= I
            push!(Ts, T)
            push!(Rs, R)
            ### print("T = $T, R = $R...")
        end
        writedlm("trans_pocket_$(string(Φ0))\_$s.dat", [Ts Rs])

        println("done after $(size(ϕ)[2]) steps.")
        writedlm("phis_pocket_im_$(string(Φ0))\_$s.dat", imag(phis))
        writedlm("phis_pocket_re_$(string(Φ0))\_$s.dat", real(phis))
    end
end

### writedlm("decay.dat", [ss Ts Rs])
end

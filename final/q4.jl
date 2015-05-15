module q4

import Evolve

steps = 8000                # number of simulation steps
samplesteps = 50            # steps between measurements
L = 2                       # width of the box
xmin = -1                   # leftmost point in the box
N = 2000                    # number of subintervals
ϵ = L/N                     # width of each subinterval
ks = [200]                  # wavenumbers
σ = 0.05                    # wavepacket width
h = 50                     # energy is 123.368 for w = 0.2
b = 125
s = 0.05                    # wall thickness
ws = [0.01:0.01:0.2]         # pocket widths to try
x0 = -0.5                   # wavepacket starting position
Ts = Float64[]
Rs = Float64[]

function pocket(L, h, b, s, w, x)
    if x < (0.5L - 0.5w - s) || x > (0.5L + 0.5w +s)
        V = 0
    elseif x < (0.5L - 0.5w) || x > (0.5L + 0.5w)
        V = h
    else
        V = b
    end
    return V
end
    
# Run the trials
for k in ks
    for w in ws
        print("Initializing for k = $k, w = $w... ")
        V(xs) = map(x -> pocket(L, h, b, s, w, x), xs)      # pocket potential
        AV = Evolve.A(V)            # evolution matrix constructor for square well

        A = AV(L, N, xmin)
        ϕ0 = zeros(Complex, N+1)
        ϕ0[int((0.5L-0.5w)N/L):int((0.5L+0.5w)N/L)] = √(10) * sin(5π*[(0.5L-0.5w):ϵ:(0.5L+0.5w)])
        x = (L/N .* [0:N] + xmin)   # x positions

        # Evolve the system and record positions every so many steps
        print("evolving... ")
        test = x -> (abs(x[2]) > 1e-2 || abs(x[end-1]) > 1e-2)
        ϕ = Evolve.evolve(A, ϕ0, steps, test)
        phis = [ϕ0 ϕ[:,samplesteps:samplesteps:end]]

        # calculate transmission and reflection probabilities
        T = mean(abs(ϕ[1:int((0.5L-0.5w-s)N/L),end]).^2)
        T += mean(abs(ϕ[int((0.5L+0.5w+s)N/L):end,end]).^2)
        R = mean(abs(ϕ[int((0.5L-0.5w)N/L):int((0.5L+0.5w)N/L),end]).^2)
        I = T + R
        T /= I
        R /= I
        push!(Ts, T)
        push!(Rs, R)
        print("T = $T, R = $R...")
        writedlm("trans_pocket_$w", [T,R,I])

        println("done after $(size(ϕ)[2]) steps. Writing to file.")
        writedlm("phis_pocket_im_$k\_$w.dat", imag(phis))
        writedlm("phis_pocket_re_$k\_$w.dat", real(phis))
    end
end

writedlm("decay.dat", [ws Ts Rs])
end

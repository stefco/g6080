module Final

import Evolve

L = 2                       # width of the box
xmin = -1                   # leftmost point in the box
N = 20000                   # number of subintervals
ϵ = L/N                     # width of each subinterval
k = 20                      # wavenumbers
σ = 0.2                     # wavepacket width
x0 = 0.0                    # wavepacket starting position

V(xs) = map(x->0.0, xs)     # zero potential energy
AV = Evolve.A(V)            # evolution matrix constructor for square well
ϕV = Evolve.Φ(k, σ, x0)     # initial state vector constructor for Gaussian packet

x = (L/N .* [0:N] + xmin)   # x positions


end

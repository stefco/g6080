include("cluster.jl")

Nx = Ny = 16
steps = 20000

Jc = 0.4406868                              # critical temp reduced coupling const
J = 1.0                                     # coupling constant in Kelvin
Tc = J/Jc                                   # critical temperature
println("Critical temperature Tc = $Tc")

Ts = [0.01:0.01:0.09, 0.1:0.1:2.3]
Bs = [0.0001, 0.001, 0.01]                     # magnetic field

Tvals = [0.001:0.001:2.269]
Mvals = M(Tvals, J)

# for plotting results
μvsT(b) = Gadfly.plot(
    Gadfly.layer(
        x=Ts, y=μs[:,b], Gadfly.Geom.point, Gadfly.Geom.line, 
        color=["Measured for B=$(Bs[b])"]),
    Gadfly.layer(
        x=Tvals, y=Mvals, Gadfly.Geom.line,
        color=["Predicted magnetization"]),
    Gadfly.Guide.xlabel("Temperature"),
    Gadfly.Guide.ylabel("Mean normalized magnetization"), 
    Gadfly.Guide.title("Magnetization vs. Temperature for B=$(Bs[b])")
)

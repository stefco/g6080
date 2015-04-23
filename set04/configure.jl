include("cluster.jl")

Nx = Ny = 16
@show steps = 2000

Jc = 0.4406868                              # critical temp reduced coupling const
J = 1.0                                     # coupling constant in Kelvin
Tc = J/Jc                                   # critical temperature
println("Critical temperature Tc = $Tc")

Ts = [2.2:0.01:2.3]
Bs = [0.0001, 0.001, 0.01]                     # magnetic field

Tvals = [0.001:0.001:2.269]
Mvals = M(Tvals, J)

print("Importing Gadfly... ")
import Gadfly
Gadfly.set_default_plot_size(20Gadfly.cm, 12Gadfly.cm)
println("done.")

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

# see magnetization distribution
μhist(mags) = Gadfly.plot(
    Gadfly.Guide.title("Magnetization Histogram for B = 0.1"),
    Gadfly.layer(
        x=abs(mags), color=["Absolute Magnetization"], 
        Gadfly.Geom.histogram(bincount=100, density=true), order=1),
    Gadfly.layer(
        x=mags, color=["Magnetization"], 
        Gadfly.Geom.histogram(bincount=100, density=true), order=2)
)

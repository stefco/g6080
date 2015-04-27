include("cluster.jl")

Nx = Ny = 128
@show steps = 40000

Jc = 0.4406868                              # critical temp reduced coupling const
J = 1.0                                     # coupling constant in Kelvin
Tc = J/Jc                                   # critical temperature in Kelvin
println("Critical temperature Tc = $Tc")

# 17 distinct temperatures
Ts = [0.5, 1.0, 1.25, 1.5:0.1:2.2, 2.22, 2.24, 2.26, 2.265, 2.27]
# magnetic fields
Bs = [0.00001, 0.0001]

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
    Gadfly.Guide.ylabel("Mean normalized absolute magnetization"), 
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

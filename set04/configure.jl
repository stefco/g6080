include("cluster.jl")

Nx = Ny = 16
@show steps = 20000

Jc = 0.4406868                              # critical temp reduced coupling const
J = 1.0                                     # coupling constant in Kelvin
Tc = J/Jc                                   # critical temperature in Kelvin
println("Critical temperature Tc = $Tc")

# Try different lattice sizes
Ns = [16, 32, 128]
# 17 distinct temperatures
Ts = [0.5, 1.0, 1.25, 1.5:0.1:2.2, 2.22, 2.24, 2.26, 2.265, 2.27, 2.28, 2.3:0.1:3.0]
# magnetic fields
Bs = [0.00001, 0.0001, 0.001, 0.01, 0.1]

Tvals = [0.001:0.001:2.269]
Mvals = M(Tvals, J)

# maximum correlation distance
mcd = 12

print("Importing Gadfly... ")
import Gadfly
Gadfly.set_default_plot_size(20Gadfly.cm, 12Gadfly.cm)
println("done.")

# for plotting results
function MeansvsT(means, b)

    layers = Array(Gadfly.Layer,0);
    for n in 1:length(Ns)
        newlayer = Gadfly.layer(
            x=Ts, y=means[:,n], Gadfly.Geom.point, Gadfly.Geom.line, 
            color=["Simulated for N=$(Ns[n])"]);
        layers = cat(1, layers, newlayer)
    end

    Gadfly.plot(
        layers,
        Gadfly.layer(
            x=Tvals, y=Mvals, Gadfly.Geom.line,
            color=["Predicted magnetization"]),
        Gadfly.layer(
            x=[Tc,Tc], y=[0,1], Gadfly.Geom.line,
            color=["Critical Temperature 2.269K"]),
        Gadfly.Guide.xticks(ticks=[0:0.5:2, 2.269, 2.5, 3]),
        Gadfly.Guide.yticks(ticks=[0:0.1:1]),
        Gadfly.Guide.xlabel("Temperature"),
        Gadfly.Guide.ylabel("Mean normalized absolute magnetization"), 
        Gadfly.Guide.title("Magnetization vs. Temperature for B=$(Bs[b]) at Different Lattice Sizes")
    )
end

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
    Gadfly.Guide.title("Magnetization vs. Temperature for B=$(Bs[b]) with $(Nx*Ny) Particles")
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

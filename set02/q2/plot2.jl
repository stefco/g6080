import Gadfly
Gadfly.set_default_plot_size(20Gadfly.cm, 12Gadfly.cm)

# make energy conservation and velocity distribution plots

temps_1069 = readdlm("temps-rho-0.75-temp-1.069.dat")
temps_1304 = readdlm("temps-rho-0.75-temp-1.304.dat")

function mdplot(values::Array, name::String, N::Integer)
    p = Gadfly.plot(
        x=[10:10:N], 
        y=temps_low, 
        Gadfly.Geom.line,
        Gadfly.Guide.xlabel("Simulation Step"), 
        Gadfly.Guide.ylabel(name),
        Gadfly.Guide.title("$name vs. Simulation Step for Liquid Argon")
    )

    # save as svg and png
    Gadfly.draw(Gadfly.SVG("$name.svg", 20Gadfly.com, 12Gadfly.cm), p)
    Gadfly.draw(Gadfly.PNG("$name.png", 12Gadfly.com, 8.cm), p)

    return p
end

function mdplottemps(temps::Array, name::String, N::Integer, T::Real)
    p = Gadfly.plot(
        layer(
            x = temps,
            Gadfly.Geom.histogram(bincount=50, density=true),
            color = ["Particle speed distribution for T = $T"]
        ),
        layer(
            x = 

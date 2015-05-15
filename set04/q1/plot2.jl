import Gadfly
Gadfly.set_default_plot_size(20Gadfly.cm, 12Gadfly.cm)

# make energy conservation and velocity distribution plots

pe = readdlm("pe.dat")
vir = readdlm("v.dat")

function mdplot(values::Array, title::String, name::String, N::Integer)
    p = Gadfly.plot(
        x=[10:10:N], 
        y=values, 
        Gadfly.Geom.line,
        Gadfly.Guide.xlabel("Simulation Step"), 
        Gadfly.Guide.ylabel(name),
        Gadfly.Guide.title("$title vs. Simulation Step for Liquid Argon")
    )

    # save as svg and png
    Gadfly.draw(Gadfly.SVG("$name.svg", 20Gadfly.cm, 12Gadfly.cm), p)
    Gadfly.draw(Gadfly.PNG("$name.png", 12Gadfly.cm, 8Gadfly.cm), p)

    return p
end

# m = 1, k = 1 in our natural units
function maxboltzbuilder(T::Real)
    return v -> √(2/π*T^3) * v.^2 .* exp(-v.^2 / (2T))
end

function mdplotspeeds(speeds::Array, name::String, T::Real)
    xs = [0:0.01:1]
    p = Gadfly.plot(
        Gadfly.layer(
            x = speeds,
            Gadfly.Geom.histogram(bincount=100, density=true),
            color = ["Measured speed histogram for T = $T"]
        ),
###         Gadfly.layer(
###             x = xs,
###             y = maxboltzbuilder(T)(xs),
###             Gadfly.Geom.line,
###             color = ["Predicted by Maxwell-Boltzmann"]
###         ),
###         color = ["Legend"],
        Gadfly.Guide.title("Temperature Distribution for T=$T vs. Boltzmann Distribution"),
        Gadfly.Guide.xlabel("Velocity"),
        Gadfly.Guide.ylabel("Probability"),
        Gadfly.Guide.xticks(ticks=[0:0.1:1])
    )

    # save as svg and png
    Gadfly.draw(Gadfly.SVG("$name.svg", 20Gadfly.cm, 12Gadfly.cm), p)
    Gadfly.draw(Gadfly.PNG("$name.png", 12Gadfly.cm, 8Gadfly.cm), p)

    return p
end


pe_plots = Gadfly.Plot[]
push!(pe_plots, mdplot(pe, "Potential Energies (T = 1.069)", "temp-1.069", 20000))

e_plots = Gadfly.Plot[]
push!(e_plots, mdplot(e_1069, "Total Energies (T = 1.069)", "temp-1.069", 20000))

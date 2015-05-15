import Gadfly
Gadfly.set_default_plot_size(20Gadfly.cm, 12Gadfly.cm)

# make energy conservation and velocity distribution plots

temps_1069 = readdlm("q2/temps-rho-0.75-temp-1.069.dat")
temps_1304 = readdlm("q2/temps-rho-0.75-temp-1.304.dat")

speeds_1069 = readdlm("q2/speeds-rho-0.75-temp-1.069.dat")
speeds_1304 = readdlm("q2/speeds-rho-0.75-temp-1.304.dat")

ke_1069 = readdlm("q2/ke-rho-0.75-temp-1.069.dat")
ke_1304 = readdlm("q2/ke-rho-0.75-temp-1.304.dat")

pe_1069 = readdlm("q2/pe-rho-0.75-temp-1.069.dat")
pe_1304 = readdlm("q2/pe-rho-0.75-temp-1.304.dat")

e_1069 = pe_1069 + ke_1069
e_1304 = pe_1304 + ke_1304

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


speed_plots = Gadfly.Plot[]
push!(speed_plots, mdplotspeeds(speeds_1069[401:end], "temp-dist-1.069", 1.069));
push!(speed_plots, mdplotspeeds(speeds_1304[401:end], "temp-dist-1.304", 1.304));

temp_plots = Gadfly.Plot[]
push!(temp_plots, mdplot(temps_1069, "Temperatures (T = 1.069)", "temp-1.069", 20000))
push!(temp_plots, mdplot(temps_1304, "Temperatures (T = 1.304)", "temp-1.304", 20000))

ke_plots = Gadfly.Plot[]
push!(ke_plots, mdplot(ke_1069, "Kinetic Energies (T = 1.069)", "temp-1.069", 20000))
push!(ke_plots, mdplot(ke_1304, "Kinetic Energies (T = 1.304)", "temp-1.304", 20000))

pe_plots = Gadfly.Plot[]
push!(pe_plots, mdplot(pe_1069, "Potential Energies (T = 1.069)", "temp-1.069", 20000))
push!(pe_plots, mdplot(pe_1304, "Potential Energies (T = 1.304)", "temp-1.304", 20000))

e_plots = Gadfly.Plot[]
push!(e_plots, mdplot(e_1069, "Total Energies (T = 1.069)", "temp-1.069", 20000))
push!(e_plots, mdplot(e_1304, "Total Energies (T = 1.304)", "temp-1.304", 20000))

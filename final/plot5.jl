using Gadfly

scatter = readdlm("scatter.dat")
scatter = [scatter; readdlm("scatter1.dat")]

w = 0.02

p = plot(
    layer(
        x = scatter[:,1].*(w/π),
        y = scatter[:,2],
        Geom.line,
        Geom.point,
        color = ["Transmission Probability"]
    ),
    layer(
        x = scatter[:,1].*(w/π),
        y = scatter[:,3],
        Geom.line,
        Geom.point,
        color = ["Reflection Probability"]
    ),
    Guide.xlabel("Multiple of Resonant Wave Number, k(w/π)"),
    Guide.ylabel("Probability"),
    Guide.title("Transmission probability vs. momentum")
)

draw(SVG("scatter.svg", 20cm, 12cm), p)

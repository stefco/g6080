using Gadfly

scatter = readdlm("scatter.dat")

p = plot(
    layer(
        x = scatter[:,1],
        y = scatter[:,2],
        Geom.line,
        Geom.point,
        color = ["Transmission Probability"]
    ),
    layer(
        x = scatter[:,1],
        y = scatter[:,3],
        Geom.line,
        Geom.point,
        color = ["Reflection Probability"]
    ),
    Guide.xlabel("Wave Number (k)"),
    Guide.ylabel("Probability"),
    Guide.title("Leakage probability through potential barrier vs. momentum of wave packet")
)

draw(SVG("scatter.svg", 20cm, 12cm), p)

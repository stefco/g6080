using Gadfly

decay = readdlm("decay.dat")

p = plot(
    layer(
        x = decay[:,1],
        y = decay[:,2],
        Geom.line,
        Geom.point,
        color = ["Escape Probability"]
    ),
    layer(
        x = decay[:,1],
        y = decay[:,3],
        Geom.line,
        Geom.point,
        color = ["Confinement Probability"]
    ),
    Guide.xlabel("Width of Pocket"),
    Guide.ylabel("Probability"),
    Guide.title("Leakage probability through potential barrier vs. width of pocket")
)

draw(SVG("decay.svg", 20cm, 12cm), p)

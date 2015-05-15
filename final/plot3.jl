using Gadfly

VTR = readdlm("VTR.dat")

p = plot(
    layer(
        x = VTR[:,1],
        y = VTR[:,2],
        Geom.line,
        Geom.point,
        color = ["Transmission Probability"]
    ),
    layer(
        x = VTR[:,1],
        y = VTR[:,3],
        Geom.line,
        Geom.point,
        color = ["Reflection Probability"]
    ),
    Guide.xlabel("Potential Energy at Step Function"),
    Guide.ylabel("Probability"),
    Guide.title("Transmission probability across step potential vs. magnitude of potential for k = 200, σ = 0.05")
)

draw(SVG("VTR.svg", 20cm, 12cm), p)

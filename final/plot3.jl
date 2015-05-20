using Gadfly

VTR = readdlm("VTR.dat")

# Central wave number for incoming wave packet
k = 200

# Theoretically predicted transmission and reflection coefficients
T = (V0) -> k^2 < 2V0 ? 0.0 : 4k*√(k^2 - 2V0) / (k+√(k^2 - 2V0))^2
R = (V0) -> k^2 < 2V0 ? 1.0 : (k-√(k^2 - 2V0))^2 / (k+√(k^2 - 2V0))^2

Vs = [5e3:100:5e4]

p = plot(
    layer(
        x = VTR[:,1],
        y = VTR[:,2],
        Geom.line,
        Geom.point,
        color = ["Transmission Probability, Experimental"]
    ),
    layer(
        x = VTR[:,1],
        y = VTR[:,3],
        Geom.line,
        Geom.point,
        color = ["Reflection Probability, Experimental"]
    ),
    layer(
        x = Vs,
        y = map(T, Vs),
        Geom.line,
        color = ["Transmission Probability, Theoretical"]
    ),
    layer(
        x = Vs,
        y = map(R, Vs),
        Geom.line,
        color = ["Reflection Probability, Theoretical"]
    ),
    Guide.xlabel("Potential Energy at Step Function"),
    Guide.ylabel("Probability"),
    Guide.title("Transmission probability across step potential vs. magnitude of potential for k = 200, σ = 0.05")
)

draw(SVG("VTR.svg", 20cm, 12cm), p)

function gadflyhistogram(bigvec, smallvec)
    using Gadfly
    plot(
        layer(
            x = bigvec,
            Geom.histogram(bincount=25),
            order=1
        ),
        layer(
            x = smallvec,
            Geom.histogram(bincount=25), 
            order=2,
            Theme(default_color=color("red"))
        ),
        Guide.xlabel("Sample mean"),
        Guide.ylabel("Frequency"),
        Guide.title("Sample Mean Histogram in Dataset 1")
    );
end

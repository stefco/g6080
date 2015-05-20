using Gadfly

Φs = [:Φground, :Φexcited]
ss = [0.001, 0.01, 0.1]                    # wall thickness

for Φ0 in Φs
    for s in ss

    # Plot the log of escape probability
    decay = readdlm("trans_pocket_$(string(Φ0))\_$s.dat")
    tstart = 5
    ts = 50*(tstart-1+[tstart:size(decay)[1]])*1e-6
    lin = linreg(ts, log(decay[tstart:end,1]))
    bestfit = x->lin[2]*x + lin[1]          # Linear Regression
    p = plot(
        layer(
            x = ts,
            y = log(decay[tstart:end,1]),
            Geom.line,
            color = ["log(T)"]
        ),
#       layer(
#           x = ts,
#           y = log(decay[tstart:end,2]),
#           Geom.line,
#           color = ["log(R)"]
#       ),
        layer(
            x = ts,
            y = bestfit(ts),
            Geom.line,
            color = ["Linear Regression of log(T), a = $(lin[2])"]
        ),
        Guide.xlabel("Number of Steps"),
        Guide.ylabel("log(Probability)"),
        Guide.title("Log of Leakage Probability vs. Time for $(string(Φ0)) initial state, s = $s")
    )
    draw(SVG("decay_$(string(Φ0))\_$s.svg", 20cm, 12cm), p)

    # Plot final vs initial configuration
    phis = readdlm("phis_pocket_re_$(string(Φ0))\_$s.dat")
    xs = 1e-3*[0:2000]-1.0
    t = 50size(phis)[2]*1e-6
    q = plot(
        layer(
            x = xs,
            y = phis[:,1],
            Geom.line,
            color = ["Wavefunction at t=0"]
        ),
        layer(
            x = xs,
            y = phis[:,end],
            Geom.line,
            color = ["Wavefunction at t=$t"]
        ),
        Guide.xlabel("x"),
        Guide.ylabel("Re(ϕ)"),
        Guide.title("ϕ initial and final starting from $(string(Φ0)) with wall thickness s = $s")
    )
    draw(SVG("wave_$(string(Φ0))\_$s.svg", 20cm, 12cm), q)

    end
end

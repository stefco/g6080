include("run2.jl")
import Gadfly
Gadfly.set_default_plot_size(20Gadfly.cm, 12Gadfly.cm)

intrvls3 = Any[
    (0,2,50),
    (36.5,0.01,37)
]

# generate data for the plots on given intervals and save it to disk
function data3(intrvl)

    for i in 1:length(intrvl)
        qs = [intrvl[i][1]:intrvl[i][2]:intrvl[i][3]]
        Ps = zeros(length(qs))
        @sync @parallel for j in 1:length(qs)
            println("Finding ψ*ψ in bottom half for q = ",qs[j])
            Ps[j] = q3P(32,qs[j],20,1e-3);
        end
        writedlm("p3probs-$(intrvl[i][1])-to-$(intrvl[i][3]).dat",Ps)
    end

end

# create plots and save them as images
function plot3(intrvl)

    plots = Array(Gadfly.Plot,0)

    @sync @parallel for i in 1:length(intrvl)
        Ps = readdlm("p3probs-$(intrvl[i][1])-to-$(intrvl[i][3]).dat")
        p = Gadfly.plot(
            x=[intrvl[i][1]:intrvl[i][2]:intrvl[i][3]],
            y=Ps,
            Gadfly.Geom.line, Gadfly.Geom.point,
            Gadfly.Guide.xlabel("Charge of Point Particle, q"),
            Gadfly.Guide.ylabel("ψ*ψ"),
            Gadfly.Guide.title("Probability of Finding Particle in Bottom Half of Region vs. Charge of Point Particle"),
            Gadfly.Guide.xticks(ticks=[intrvl[i][1]:5intrvl[i][2]:intrvl[i][3]])
        )
        # save as SVG
        Gadfly.draw(Gadfly.SVG("p3probs-$(intrvl[i][1])-to-$(intrvl[i][3]).svg", 20Gadfly.cm, 12Gadfly.cm), p)
        push!(plots,p)

    end

    return plots

end

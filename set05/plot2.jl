include("run1.jl")
import Gadfly
# Gadfly.set_default_plot_size(20Gadfly.cm, 12Gadfly.cm)

intrvl = {
    (1,1,600),
    (153,0.01,154),
    (272,0.01,273),
    (368,0.01,369)
}

# generate data for the plots on given intervals
function data2(intrvl)

    datasets = Array(Array{Float64},length(intrvl))

    @sync @parallel for i in 1:length(intrvl)
        eigval, Δe = q2(32, [intrvl[i][1]:intrvl[i][2]:intrvl[i][3]]);
        datasets[i] = eigval
        writedlm("p2eigvals-$(intrvl[i][1])-to-$(intrvl[i][3]).dat",eigval)
    end

    return datasets

end

# create plots and save them as images
function plot2(intrvl)

    @sync @parallel for i in 1:length(intrvl)
        eigval = readdlm("p2eigvals-$(intrvl[i][1])-to-$(intrvl[i][3]).dat")
        p = Gadfly.plot(
            x=[intrvl[i][1]:intrvl[i][2]:intrvl[i][3]],
            y=eigval,
            Gadfly.Geom.line, Gadfly.Geom.point,
            Gadfly.Guide.xlabel("Eigenvalue Estimate"), 
            Gadfly.Guide.ylabel("Magnitude of ψ*ψ"),
            Gadfly.Guide.title("Magnitude of wavefunction squared vs. Estimated Eigenvalue"),
            Gadfly.Guide.xticks(ticks=[intrvl[i][1]:10intrvl[i][2]:intrvl[i][3]])
        )
        # save as SVG
        Gadfly.draw(Gadfly.SVG("p2eigvals-$(intrvl[i][1])-to-$(intrvl[i][3]).svg", 20Gadfly.cm, 12Gadfly.cm), p)

    end

end

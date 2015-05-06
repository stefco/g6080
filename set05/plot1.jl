import Gadfly
include("run1.jl")
Gadfly.set_default_plot_size(20Gadfly.cm, 12Gadfly.cm)

pts1 = [(2,2), (3,2), (4,4), (7,7), (7,5), (7,1)]

function data1(pts::Array)
    # Set of points, in parallel
    vals = Array(Float64, 8, length(pts))

    # calculate the potential for each of several grid sizes
    @sync @parallel for M in 1:8
        ϕ = q1(8M)[1]               # only grab the potential
        for i in 1:length(pts)
            x = M*pts[i][1] + 1
            y = M*pts[i][2] + 1
            vals[M,i] = ϕ[x,y]
        end

    end

    # must write to file
    for i in 1:length(pts)
        writedlm("p1phi-$(pts[i][1]/8)-$(pts[i][2]/8).dat", vals[:,i])
    end

    return vals

end

function plot1(pts)

    # Initialize containers with space for n layers for each point
    plots = Array(Gadfly.Plot, length(pts))

    # For each point, put all the layers in a Gadfly plot with appropriate legend
    @sync @parallel for i in 1:length(pts)
        ϕs = readdlm("p1phi-$(pts[i][1]/8)-$(pts[i][2]/8).dat")
        p = Gadfly.plot(
            x=(8*[1:8]),
            y=ϕs,
            Gadfly.Geom.line, Gadfly.Geom.point,
            Gadfly.Guide.xlabel("Horizontal Divisions, N"), 
            Gadfly.Guide.ylabel("ϕ"),
            Gadfly.Guide.title("Magnitude of ϕ vs. Horizontal Grid Divisions at x = $(pts[i][1]/8), y = $(pts[i][2]/8)"),
            Gadfly.Guide.xticks(ticks=(8*[1:8]))
        )
        # save as SVG
        Gadfly.draw(Gadfly.SVG("p1phi-$(pts[i][1]/8)-$(pts[i][2]/8).svg", 20Gadfly.cm, 12Gadfly.cm), p)
        plots[i] = p

    end

    return plots
end

# Include everything
include("md.jl")

# File to load and save to
fname = "../../swp/s_auto.mdv"

# Load trial s_auto
s = loadmdtrial(fname)

# Run for another 10000 steps
@time for i in 1:20
    println("starting another round")
    verlet!(s, 500)
    save(fname, s)
    println("saved round $i.")
end

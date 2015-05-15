## QUESTION 2

include("md.jl");

# Pick interaction distance
δ = 0.1

# Create a new simulation object
qtherm = MolecularDynamicsTrial(
    1000,               # number of bodies
    0.75,               # particle density
    1.069,              # desired temperature; also use 1.304
    0.032,              # time step dx, in units of τ
    2000,              # number of iterations
    3.405e-10,          # length units, from Verlet, in m
    119.8,              # energy scale for potential, in Kelvin
    39.948 / 6.0221e26  # argon atomic mass, in kg
);      

# Place the particles
mdplace!(qtherm);

# Run for 2000 steps
metropolis!(qtherm, 1999, δ)

# Save the trial
print("Done thermalizing. Saving new starting position...")
save("../../swp/1069_met_therm.mdv", qtherm)

# Take the last step as the start step
qstart = subtrial(qtherm, 2000, 2000)
save("../../swp/1069_met_start.mdv", qstart)
println("done.")

## QUESTION 2

include("md.jl");

# Pick interaction distance
δ = 0.1

# Create a new simulation object
q = MolecularDynamicsTrial(
    1000,               # number of bodies
    0.75,               # particle density
    1.069,              # desired temperature; also use 1.304
    0.032,              # time step dx, in units of τ
    10000,              # number of iterations
    3.405e-10,          # length units, from Verlet, in m
    119.8,              # energy scale for potential, in Kelvin
    39.948 / 6.0221e26  # argon atomic mass, in kg
);      

# Place the particles
mdplace!(q);

# Run for 199 steps
metropolis!(q, 499, δ)

# Save the trial
save("set04/q1/q.mdv", q)

# WILL NEED TO ADD A THERMALIZE STEP BACK IN AT SOME POINT

# Run for another 8000 steps, saving all the way
for i in 1:19
    metropolis!(q, 500, δ)
    save("set04/q1/q.mdv", q)
end

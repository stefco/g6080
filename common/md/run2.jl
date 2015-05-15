## QUESTION 2

include("md.jl");

fname = "../../swp/1304_auto.mdv"

# Create a new simulation object
q = MolecularDynamicsTrial(
    1000,               # number of bodies
    0.75,               # particle density
    1.304,              # desired temperature; also use 1.304
    0.032,              # time step dx, in units of τ
    20000,              # number of iterations
    3.405e-10,          # length units, from Verlet, in m
    119.8,              # energy scale for potential, in Kelvin
    39.948 / 6.0221e26  # argon atomic mass, in kg
);      

# Place the particles
mdplace!(q);

# Run for 199 steps to thermalize
verlet!(q, 199)

# Save the trial
save(fname, q)

# Relax
relax!(q)

# Save the relaxed trial
save(fname, q)

# Thermalize and relax 19 more times
thermalize(fname, q, 200, 19)

# Run for another 16000 steps
for i in 1:32
    verlet!(q, 500)
    save(fname, q)
end

# Extract information
extract("../../swp/set02", q, 10)

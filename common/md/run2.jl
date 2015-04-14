## QUESTION 2

include("md.jl");

# Create a new simulation object
q = MolecularDynamicsTrial(
    1000,               # number of bodies
    0.75,               # particle density
    1.304,              # desired temperature; also use 1.304
    0.032,              # time step dx, in units of τ
    10000,              # number of iterations
    3.405e-10,          # length units, from Verlet, in m
    119.8,              # energy scale for potential, in Kelvin
    39.948 / 6.0221e26  # argon atomic mass, in kg
    );      

# Create a new index object
index = MolecularDynamicsIndex( 3, 10 )

# Place the particles
mdplace!(q, index);

# Run for 199 steps to thermalize
verlet!(q, 199, index)

# Save the trial
save("set02/q2/q200", q)

# Relax
relax!(q)

# Save the relaxed trial
save("set02/q2/q200r", q)

# Thermalize
thermalize("set02/q2/q", q, 50, 36, index)

# Run for another 8000 steps
verlet!(q, index)

# Save
save("set02/q2/qfinal",q)

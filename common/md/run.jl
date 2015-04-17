## QUESTION 2

include("md.jl");

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

# Run for 199 steps to thermalize
verlet!(q, 199)

# Save the trial
save("set02/q2/q_auto", q)

# Relax
relax!(q)

# Save the relaxed trial
save("set02/q2/q_auto", q)

# Thermalize
thermalize("set02/q2/q_auto", q, 50, 36)

# Run for another 8000 steps
for i in 1:16
    verlet!(q, 500)
    save("set02/q2/q_auto", q)
end

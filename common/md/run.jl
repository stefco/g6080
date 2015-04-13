## QUESTION 2

include("MolecularDynamicsTrial.jl");
include("mdplace.jl");

# Create a new simulation object
r = MolecularDynamicsTrial(
    1000,               # number of bodies
    0.75,               # particle density
    1.069,              # desired temperature; also use 1.304
    0.032,              # time step dx, in units of τ
    1e5,                # number of iterations
    3.405e-10,          # length units, from Verlet, in m
    119.8,              # energy scale for potential, in Kelvin
    39.948 / 6.0221e26, # argon atomic mass, in kg
    );      

# Place the particles
mdplace!(r);



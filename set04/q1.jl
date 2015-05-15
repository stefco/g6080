module q1

# function for calculating temperature and pressure
function tandp(ke, pe, vir, N)
    # Take the time average of the virial, ⟨∑∑r⃗∂V/∂r⟩
    virialavg = means(vir)

    # Find the temperatures
    T = 2ke/3N

    # Set Boltsmann constant
    kB = 1.3806e-23

    kBTs = kB*T

    p = 1 - virialavg./(6N*kBTs)

    return T, p
end

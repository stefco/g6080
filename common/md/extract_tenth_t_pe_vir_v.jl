# take a given molecular dynamics trial and save every tenth PE, KE, T, virial, and 
# speed distribution to a file whose name is based on ρ and Tdesired.
function extract(r::MolecularDynamicsTrial, N::Integer)
    rho = r.ρ
    Td  = r.Td

    speeds = √sum(r.v[:,:,N:N:end] .* r.v[:,:,N:N:end], 1)
    speeds = reshape(speeds, size(speeds)[2], size(speeds)[3])
    writedlm("speeds-rho-$rho-temp-$Td.dat", speeds)

    temps = r.T[N:N:end]
    writedlm("temps-rho-$rho-temp-$Td.dat", temps)

    ke = r.ket[N:N:end]
    writedlm("ke-rho-$rho-temp-$Td.dat", ke)

    pe = r.pet[N:N:end]
    writedlm("pe-rho-$rho-temp-$Td.dat", pe)

    vir = r.vir[N:N:end]
    writedlm("vir-rho-$rho-temp-$Td.dat", vir)
end

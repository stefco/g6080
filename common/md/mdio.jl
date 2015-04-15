# Methods for manipulating MolecularDynamicsTrial objects

include("MolecularDynamicsTrial.jl")
import JSON
import HDF5

# Write the trial to JSON
function writejson( filename::String, r::MolecularDynamicsTrial )
    isfile( filename ) && rm( filename )
    file = open( filename, "w" )
    write( file, JSON.json( r, 4 ) )
    close( file )
end

# Save the trial using HDF5
function save( filename::String, r::MolecularDynamicsTrial )
    print("Saving trial as ",filename,"...")
    fileid = HDF5.h5open( filename, "w" )               # open, overwrite
    for field in names( MolecularDynamicsTrial )
        fileid[ string( field ) ] = r.( field )
    end
    close( fileid )
    println(" done.")
end

# Load the trial using HDF5
function load!( filename::String, r::MolecularDynamicsTrial )
    print("Loading trial from ",filename,"...")
    fileid = HDF5.h5open( filename, "r" )               # open, read-only
    for field in names( fileid )
        r.( symbol( field ) ) = HDF5.read( fileid, field )
    end
    close( fileid )
    println(" done.")
end

# Initialize and load a trial using HDF5
function loadmdtrial( filename::String )
    r = MolecularDynamicsTrial()
    load!( filename, r )
    return r
end

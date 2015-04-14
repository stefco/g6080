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
    fileid = HDF5.h5open( filename, "w" )               # open, overwrite
    for field in names( MolecularDynamicsTrial )        # save all fields
        fileid[ string( field ) ] = r.( field )
    end
    close( fileid )                                     # close and save file
end

# Load the trial using HDF5
function load!( filename::String, r::MolecularDynamicsTrial )
    fileid = HDF5.h5open( filename, "r" )               # open, read-only
    for field in names( fileid )                        # populate fields
        r.( symbol( field ) ) = HDF5.read( fileid, field )
    end
    close( fileid )                                     # close file
end

# Initialize and load a trial using HDF5
function loadmdtrial( filename::String )
    r = MolecularDynamicsTrial()                        # create new trial
    load!( filename, r )                                # populate from save
    return r
end

using LaMEM, GeophysicalModelGenerator

# structure definition
struct markerGrid
    Xpart::Array{Float64,3}
    Ypart::Array{Float64,3}
    Zpart::Array{Float64,3}
    Temp::Array{Float64,3}
    Phase::Array{Int64,3}
end

markerGrid(X, Y, Z) = markerGrid(X, Y, Z, zeros(size(X)), zeros(Int64, size(X)))


# function definition
function createMarkerGrid(filename::String, numCores::Int64)
    # create partitioning file from .dat file
    partName = run_lamem_save_grid(filename, numCores)

    # read .dat file to make grid
    Grid     = ReadLaMEM_InputFile(filename)

    # create marker Grid
    A        = markerGrid(Grid.X, Grid.Y, Grid.Z)

    return A, partName
end
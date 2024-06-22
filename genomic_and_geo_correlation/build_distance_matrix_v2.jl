"""
    Build geographical distance matrices
    This version is (much) slower but works offline with downloaded maps.
    It does not require access to the OSRM api.
    This script serves as an example of how one could proceed with download maps,
    but will thus not be tested/maintained!
"""

# Author: Neil Scheidwasser (neil.clow@sund.ku.dk)

using ArgParse
using Base.Threads
using CSV
using Combinatorics
using DataFrames
using Dates
using JLD2
using OpenStreetMapX
using ProgressMeter
using Random

# Aliases
ox = OpenStreetMapX

# Paths
GEO_PATH = "data/synthetic/geo/"

MAP_PATH = "data/maps"

# Functions
"""
	parse_commandline()

Parse command line arguments.

# Returns
- `Dict{String, Any}`: A dictionary containing parsed command line arguments.
"""
function parse_commandline()
    args = ArgParseSettings()

    @add_arg_table args begin
        "--location", "-L"
        arg_type = String
        help = "Location from which data is taken (nothing: all of DK)"
    end

    return parse_args(args)
end


"""
	load_data(location)

Load local, regional or national metadata.

# Arguments
- `location::String`: The name of the location to analyse.

# Returns
- `df::DataFrame`: A DataFrame containing the loaded data
- `lon::Vector{Float64}`: Vector of longitudes.
- `lat::Vector{Float64}`: Vector of latitudes.
- `map_fname::String`: Path to a map corresponding to the locaion of intereste.
"""
function load_data(location)
    if isnothing(location)
        location = "national"
    end

    # Find corresponding map
    if location in ["Hovedstaden", "Midtjylland", "Nordjylland", "Sjælland", "Syddanmark"]
        map_fname = "$MAP_PATH/$(lowercase(location)).osm.pbf"
    elseif location == "Copenhagen"
        map_fame = "$MAP_PATH/Hovedstaden.osm.pbf"
    elseif location == "Aarhus"
        map_fame = "$MAP_PATH/Midtjylland.osm.pbf"
    elseif location == "Odense"
        map_fame = "$MAP_PATH/Syddanmark.osm.pbf"
    else
        map_fname = "$MAP_PATH/denmark-latest.osm.pbf"
    end

    # Load metadata
    metadata_fname = "$GEO_PATH/$location.csv"

    df = DataFrame(CSV.File(metadata_fname))

    println("Number of samples: $(size(df, 1))")

    lon, lat = df[!, "longitude"], df[!, "latitude"]

    return df, lon, lat, map_fname
end

"""
	build_map(map_fname)

Build a network map from a filename.

# Arguments
- `map_fname::String`: The filename of the map data.

# Returns
- `OpenStreetMapX::MapData`: Map data.
"""
function build_map(map_fname)
    map_data = ox.get_map_data(map_fname)

    return map_data
end

"""
	get_nearest_nodes(map_data, sub_df, lon, lat, n_samples, save=true)

Get the nearest nodes to given coordinates.

# Arguments
- `map_data::Any`: Map data.
- `df::DataFrame`: Subset of data containing coordinates.
- `lon::Vector{Float64}`: Vector of longitudes.
- `lat::Vector{Float64}`: Vector of latitudes.
- `n_samples::Int`: Number of samples.
- `save::Bool`: Whether to save the result. Default is `true`.

# Returns
- `Vector{Int}`: Vector of nearest nodes.
"""
function get_nearest_nodes(map_data, df, lon, lat, n_samples, save=true)
    nearest_nodes = Vector{Int}(undef, length(lon))

    p = Progress(n_samples, showspeed=true)

    Threads.@threads for (i, (lat_i, lon_i)) in collect(enumerate(zip(lat, lon)))
        nearest_nodes[i] = ox.point_to_nodes((lat_i, lon_i), map_data)

        next!(p)
    end

    df[!, :nearest_node] = nearest_nodes

    sort!(df, [:nearest_node])

    if save
        CSV.write("data/geo_metadata-n$n_samples.csv")
    end

    return df.nearest_node
end

"""
	build_distance_matrix(map_data, nearest_nodes, n_samples)

Build a distance matrix from map data and nearest nodes.

# Arguments
- `map_data::Any`: Map data.
- `nearest_nodes::Vector{Int}`: Vector of nearest nodes.
- `n_samples::Int`: Number of samples.

# Returns
- `Vector{Float64}`: Route distances.
- `Vector{Float64}`: Route times.
"""
function build_distance_matrix(map_data, nearest_nodes, n_samples)
    route_distances = zeros(n_samples, n_samples)
    route_times = zeros(n_samples, n_samples)

    n_steps = Int(0.5 * n_samples * (n_samples - 1))

    p = Progress(n_steps, showspeed=true)

    Threads.@threads for (i, j) in collect(combinations(1:n_samples, 2))
        next!(p)

        _, route_distance, route_time = ox.shortest_route(map_data, nearest_nodes[i], nearest_nodes[j])
        route_distances[i, j] = route_distance
        route_times[i, j] = route_time
    end

    return route_distances, route_times

end


function main()
    args = parse_commandline()

    println("Configuration:")
    for (k, v) in args
        println("\t$k: $v")
    end
    println("\tThreads: $(Threads.nthreads())")

    location = args["location"]

    println("Loading data...")
    df, lon, lat, map_fname = load_data(location)
    n_samples = nrow(df)

    println("Loading map...")
    map_data = build_map(map_fname)

    # Prepare outputs
    distance_folder = "data/distance_matrices/distances"
    time_folder = "data/distance_matrices/time"

    mkpath(distance_folder)
    mkpath(time_folder)

    distances_output = "$distance_folder/route_distances-l$location-n$(length(lon)).jld2"
    times_output = "$time_folder/route_times-l$location-n$(length(lon)).jld2"

    println("Finding nearest nodes...")
    nearest_nodes = get_nearest_nodes(map_data, df, lon, lat, n_samples)

    println("Building distance matrix...")
    route_distances, route_times = build_distance_matrix(map_data, nearest_nodes, n_samples)

    println("Saving distance matrices...")
    save(
        distances_output,
        "route_distances",
        route_distances,
    )
    save(
        times_output,
        "route_times",
        route_times,
    )
    println("Done")
end

# Call the main function if this script is executed directly
if abspath(PROGRAM_FILE) == @__FILE__
    main()
end

"""Build geographical distance matrices"""

# Author: Neil Scheidwasser (neil.clow@sund.ku.dk)

using ArgParse
using Base.Threads
using CSV
using Combinatorics
using DataFrames
using HTTP
using JLD2
using JSON
using ProgressMeter

GEO_PATH = "data/synthetic/geo/"


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
        "--fill", "-f"
        action = :store_true
        help = "If true, only fill NaNs and zeros in already saved matrices"
    end

    return parse_args(args)
end

"""
	load_data(location)

Load local, regional or national metadata.

# Arguments
- `location::String`: The name of the location to analyse.

# Returns
- `lon::Vector{Float64}`: Vector of longitudes.
- `lat::Vector{Float64}`: Vector of latitudes.
"""
function load_data(location)
    if isnothing(location)
        location = "national"
    end

    metadata_fname = "$GEO_PATH/$location.csv"

    df = DataFrame(CSV.File(metadata_fname))

    println("Number of samples: $(size(df, 1))")

    lon, lat = df[!, "lon"], df[!, "lat"]

    return lon, lat
end


"""
	get_distance(lon, lat, i, j)

Compute the distance between two points given their longitude and latitude coordinates.

# Arguments
- `lon::Vector{Float64}`: Vector of longitudes.
- `lat::Vector{Float64}`: Vector of latitudes.
- `i::Int64`: The index of the first point.
- `j::Int64`: The index of the second point.

# Returns
- `Float64, Float64`: route distance and route time between first and second point.
"""
function get_distance(lon, lat, i, j)
    try
        response = HTTP.get("http://localhost:5000/route/v1/driving/$(lon[i]),$(lat[i]);$(lon[j]),$(lat[j])", readtimeout=1)
        if response.status == 200
            stats = JSON.parse(String(response.body))["routes"][1]
            return stats["distance"], stats["duration"]
        else
            return NaN, NaN
        end
    catch ex
        if isa(ex, Union{HTTP.Exceptions.TimeoutError,HTTP.Exceptions.ConnectError,HTTP.Exceptions.StatusError})
            # Timeout and connect errors: bugs from osrm api or apptainer
            # Status error: probably for impossible distances involving ferry travel (e.g., Bornholm)
            return NaN, NaN
        else
            rethrow(ex)
        end
    end
end

"""
	load_data(lon, lat)

Build a distance a matrix from longitude-latitude vectors.

Uses the OSRM API.

# Arguments
- `lon::Vector{Float64}`: Vector of longitudes (size: [n_samples, 1]).
- `lat::Vector{Float64}`: Vector of latitudes (size: [n_samples, 1]).

# Returns
- `route_distances::Vector{Float64, Float64}`: driving distance matrix (size: [n_samples, n_samples]).
- `route_times::Vector{Float64, Float64}`: driving time matrix (size: [n_samples, n_samples]).
"""
function build_distance_matrix(lon, lat)
    n_samples = length(lon)
    route_distances = zeros(n_samples, n_samples)
    route_times = zeros(n_samples, n_samples)
    n_steps = Int(0.5 * n_samples * (n_samples - 1))

    println("Number of steps: $n_steps")

    p = Progress(n_steps, showspeed=true)
    Threads.@threads for (i, j) in collect(combinations(1:n_samples, 2))
        next!(p)

        dist_ij, time_ij = get_distance(lon, lat, i, j)

        route_distances[i, j] = dist_ij
        route_times[i, j] = time_ij
    end

    return route_distances, route_times
end

"""
	fill_distance_matrix!(route_distances, route_times, lon, lat)

Fill NaN and zero values in route_distances and route_times

Uses the OSRM API.

NaNs can occur during build_distance_matrix for various connection
or status errors while caling the OSRM API, so this is a "second pass"
to make sure the matrices are correct.

# Arguments
- `route_distances::Vector{Float64, Float64}`: driving distance matrix (size: [n_samples, n_samples]).
- `route_times::Vector{Float64, Float64}`: driving time matrix (size: [n_samples, n_samples]).
- `lon::Vector{Float64}`: Vector of longitudes (size: [n_samples, 1]).
- `lat::Vector{Float64}`: Vector of latitudes (size: [n_samples, 1]).
"""
function fill_distance_matrix!(route_distances, route_times, lon, lat)
    nan_idxs = findall(isnan, route_distances)

    triu_zero_idxs = filter(y -> y[1] < y[2], findall(x -> x == 0, route_distances))

    idxs_to_fill = [nan_idxs; triu_zero_idxs]

    n_steps = length(idxs_to_fill)

    println("Number of steps: $n_steps")

    p = Progress(n_steps, showspeed=true)
    Threads.@threads for elm in idxs_to_fill
        next!(p)

        dist_ij, time_ij = get_distance(lon, lat, elm[1], elm[2])

        route_distances[elm] = dist_ij
        route_times[elm] = time_ij
    end
end

"""
	main()

Main script:
 - Parse arguments
 - Load metadata
 - Build distance matrices or fill/clean pre-existing ones
"""
function main()
    args = parse_commandline()

    println("Configuration:")
    for (k, v) in args
        println("\t$k: $v")
    end
    println("\tThreads: $(Threads.nthreads())")

    location = args["location"]
    do_fill = args["fill"]

    println("Loading data...")
    lon, lat = load_data(location)

    distance_folder = "data/distance_matrices/distances"
    time_folder = "data/distance_matrices/time"

    mkpath(distance_folder)
    mkpath(time_folder)

    distances_output = "$distance_folder/route_distances-l$location-n$(length(lon)).jld2"
    times_output = "$time_folder/route_times-l$location-n$(length(lon)).jld2"

    if do_fill
        println("Filling NaNs and potential erroneous zeros...")
        route_distances, route_times = load(distances_output, "route_distances"), load(times_output, "route_times")
        fill_distance_matrix!(route_distances, route_times, lon, lat)
    else
        println("Building distance matrix...")
        route_distances, route_times = build_distance_matrix(lon, lat)
    end

    nan_count = count(isnan, route_distances)
    println("NaN values: $nan_count")

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

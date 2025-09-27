function DrWatson.datadir(args...; on_cluster::Bool=false)
    if !on_cluster
        joinpath(projectdir(), "data", args...)
    else # if running on cluster data should be saved to scratch
        joinpath(ENV["SCRATCH"], projectname(), args...)
    end
end

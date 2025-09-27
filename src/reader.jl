export collect_datasets

"""
    collect_datasets(folder; rinclude, rexclude, kwargs...)
Read CSV files in folder and load them as a single dataframe.
Parameter values in the filenames are included as dataframe columns.
The `rinclude` and `rexclude` keywords can be used to filter filenames:
only datasets that match all `rinclude` regexes and none of the `rexclude`
are loaded.

Additional kwargs can be used as a more convenient way to filter
datasets based on parameter values.
Any kwarg of the form `key=val` is processed as a regex
so that only datasets where the parameter `key` has value `val`
will be loaded.
"""
function collect_datasets(folder;
    rinclude::AbstractVector{Regex}=[r""],
    rexclude::AbstractVector{Regex}=[r"^\b$"], 
    kwargs...
)
    valid_filetypes = [".csv"]
    load_function = (filename) -> CSV.read(joinpath(folder, filename), DataFrame)
    filelist = readdir(folder)
    filter!(filename -> is_valid_file(filename, valid_filetypes), filelist)
    if !isempty(kwargs)
        for (k,v) in zip(keys(kwargs), values(kwargs))
            push!(rinclude, Regex("$(k)=$(v)[^a-zA-Z0-9]"))
        end
    end
    if (rinclude == [r""] && rexclude == [r"^\b$"]) == false
        idx_filt = Int[]
        for i in eachindex(filelist)
            file = filelist[i]
            include_bool = all(!isnothing(match(rgx, file)) for rgx in rinclude)
            exclude_bool = all(!isnothing(match(rgx, file)) for rgx in rexclude)
            if !include_bool || exclude_bool
                push!(idx_filt, i)
            end
        end
        deleteat!(filelist, idx_filt)
    end
    datasets = load_function.(filelist)
    for i in eachindex(datasets)
        config = DataFrame(parse_savename(filelist[i])[2])
        datasets[i] = hcat(datasets[i], repeat(config, nrow(datasets[i])))
    end
    vcat(datasets...)
end

is_valid_file(file, valid_filetypes) = 
    any(endswith(file, v) for v in valid_filetypes)

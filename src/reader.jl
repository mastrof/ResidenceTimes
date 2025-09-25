export collect_datasets

function collect_datasets(folder; rinclude=[r""], rexclude=[r"^\b$"])
    valid_filetypes = [".csv"]
    load_function = (filename) -> CSV.read(joinpath(folder, filename), DataFrame)
    filelist = readdir(folder)
    filter!(filename -> is_valid_file(filename, valid_filetypes), filelist)
    if (rinclude == [r""] && rexclude == [r"^\b$"]) == false
        idx_filt = Int[]
        for i in eachindex(filelist)
            file = filelist[i]
            include_bool = any(!isnothing(match(rgx, file)) for rgx in rinclude)
            exclude_bool = any(!isnothing(match(rgx, file)) for rgx in rexclude)
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

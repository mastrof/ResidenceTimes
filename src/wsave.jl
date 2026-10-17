function DrWatson._wsave(
    filename,
    data::AbstractDataFrame,
    args...;
    kwargs...
)
    if endswith(filename, ".arrow")
        Arrow.write(filename, data; kwargs...)
    else
        CSV.write(filename, data; kwargs...)
    end
end

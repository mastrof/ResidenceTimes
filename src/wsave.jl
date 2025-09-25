function DrWatson._wsave(
    filename,
    data::AbstractDataFrame,
    args...;
    kwargs...
)
    CSV.write(filename, data; kwargs...)
end

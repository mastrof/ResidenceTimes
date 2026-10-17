using Distributed
@everywhere using DrWatson
@everywhere @quickactivate "ResidenceTimes"
@everywhere begin
    using ResidenceTimes
    using CSV
    using DataFrames
    using Agents
    using MicrobeAgents
end

# replay the coarse pass from its end-of-equilibration checkpoint
# and track the selected bacteria (see community_best-worst-performers.jl)
# at every timestep; positions are unfolded on the fly
@everywhere function run_replay(selection_file; on_cluster)
    ids = CSV.read(selection_file, DataFrame).id
    statefile(tag) = datadir(
        "comm_7classes_states",
        replace(basename(selection_file), r"^sel" => "state-$tag", r"\.csv$" => ".jls");
        on_cluster
    )
    model = load_checkpoint(statefile("eq"))
    final = load_checkpoint(statefile("end"))
    nsteps = abmtime(final) - abmtime(model)
    k = length(ids)
    L = spacesize(model)
    time = Vector{Int32}(undef, nsteps*k)
    id = Vector{Int32}(undef, nsteps*k)
    x = Vector{Float64}(undef, nsteps*k)
    y = Vector{Float64}(undef, nsteps*k)
    z = Vector{Float64}(undef, nsteps*k)
    c = Vector{Float64}(undef, nsteps*k)
    # excess concentration (C - Cb), same quantity that builds up the exposure
    measured = model.neighborlist.measurements.c
    previous = [position(model[i]) for i in ids]
    unfolded = copy(previous)
    for s in 1:nsteps
        step!(model, 1)
        for (j, i) in enumerate(ids)
            n = (s-1)*k + j
            pos = position(model[i])
            Δ = pos - previous[j]
            unfolded[j] += Δ - L .* round.(Δ ./ L)
            previous[j] = pos
            time[n] = abmtime(model)
            id[n] = i
            x[n], y[n], z[n] = unfolded[j]
            c[n] = measured[i]
        end
    end
    # the replay must land exactly on the coarse pass final state
    identical = isidentical(model, final)
    identical || @warn "replay diverged from the coarse pass" selection_file
    DataFrame(; time, id, x, y, z, c), identical
end

on_cluster = haskey(ENV, "SCRATCH")
selections = filter(
    endswith(".csv"),
    readdir(datadir("comm_7classes_selection"; on_cluster); join=true)
)
pmap(selections) do selection_file
    @show selection_file
    fileout = datadir(
        "comm_7classes_replay",
        replace(basename(selection_file), r"^sel" => "tracks", r"\.csv$" => ".arrow");
        on_cluster
    )
    isfile(fileout) && return
    data, identical = run_replay(selection_file; on_cluster)
    wsave(fileout, data; metadata=["replay_identical" => string(identical)])
end

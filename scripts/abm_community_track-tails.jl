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

@everywhere begin
    # used for collecting data
    _x(a::AbstractMicrobe) = position(a)[1]
    _y(a::AbstractMicrobe) = position(a)[2]
    _z(a::AbstractMicrobe) = position(a)[3]
end

parameters = Dict(
    :dt => 0.1, # s
    :mot => ["RRF"]#, "RT"],
    :U => [10, 40], # μm/s
    :λ => 2.2, # 1/s
    :PER => 0.5,
    :Cb => 0.0, # μM
    :α => 1.0,
    :Aphy => Int.([3e4, 3e5]), # cells/mL
    :L => [ # μm
        @onlyif(:Aphy == 3e4, 2150),
        @onlyif(:Aphy == 3e5, 1000),
    ],
)
dicts = dict_list(parameters)

@everywhere function run_abm(config::Dict)
    @unpack mot, dt, U, λ, L, PER, Cb, α, Aphy = config
    # read ids of 1% best and worst performers
    # run simulation with same rng seed
    # track full trajectory of only those bacteria
    ids_to_track = CSV.read(datadir(
        "comm_7classes_ids",
        savename("ids1-99", config, "csv")
    ), DataFrame).id .|> Int
    γ = 150
    model = setup_abm_community(; mot, dt, U, λ, L, α, Aphy, PER, Cb, γ, n=100_000)
    simtime = 60 # minutes
    nsteps = round(Int, simtime * 60 / dt)
    # data collection routines
    c(m::ABM, ids) = [concentration(m)(a, m) - Cb for a in allagents(m) if a.id ∈ ids]
    x(m::ABM, ids) = [_x(a) for a in allagents(m) if a.id ∈ ids]
    y(m::ABM, ids) = [_y(a) for a in allagents(m) if a.id ∈ ids]
    z(m::ABM, ids) = [_z(a) for a in allagents(m) if a.id ∈ ids]
    c(m::ABM) = c(m, ids_to_track)
    x(m::ABM) = x(m, ids_to_track)
    y(m::ABM) = y(m, ids_to_track)
    z(m::ABM) = z(m, ids_to_track)
    mdata = [x, y, z, c]
    # do not collect during first 15 minutes equilibration
    # then collect every timestep (0.1 s)
    eqtime = round(Int, 15*60 / dt)
    when_model(model, t) = t >= eqtime
    _, mdf = run!(model, nsteps; mdata, when_model, init=false)
    # rename because local functions don't respect names
    colnames = names(mdf)
    rename!(mdf,
        Symbol(colnames[2]) => :x,
        Symbol(colnames[3]) => :y,
        Symbol(colnames[4]) => :z,
        Symbol(colnames[5]) => :c,
    )
    # convert the mdf back into adf format
    adf = DataFrame(;
        time=Int64[],
        id=Int64[],
        position=SVector{3,Float64}[],
        c=Float64[],
    )
    for row in eachrow(mdf)
        time = fill(row.time, length(ids_to_track))
        id = ids_to_track
        x = row.x
        y = row.y
        z = row.z
        position = SVector{3,Float64}.(zip(x,y,z))
        c = row.c
        append!(adf, DataFrame(; time, id, position, c))
    end
    # unfold periodic positions
    Analysis.unfold!(adf, model)
    # expand unfolded position → x, y, z
    transform!(adf, :position_unfold => identity => [:x, :y, :z])
    select!(adf, Not(:position, :position_unfold))
    adf
end

pmap(dicts) do config
    @show config
    data = run_abm(config)
    on_cluster = haskey(ENV, "SCRATCH")
    fileout = datadir(
        "comm_7classes_tracks",
        savename("comm", config, "csv");
        on_cluster
    )
    wsave(fileout, data)
end

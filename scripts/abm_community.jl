using Distributed
@everywhere using DrWatson
@everywhere @quickactivate "ResidenceTimes"
@everywhere begin
    using ResidenceTimes
    using DataFrames
    using MicrobeAgents
end

@everywhere begin
    # used for collecting data
    x(a::AbstractMicrobe) = position(a)[1]
    y(a::AbstractMicrobe) = position(a)[2]
    z(a::AbstractMicrobe) = position(a)[3]
end

parameters = Dict(
    :dt => 0.1, # s
    :mot => "RRF",
    :U => [10, 20, 40], # μm/s
    :λ => 2.2, # 1/s
    :PER => [0.3, 0.5],
    :Cb => 0.0, # μM
    :α => [0.75, 1.0],
    :Aphy => Int.([3e4, 8e4, 3e5]), # cells/mL
    :L => [ # μm
        @onlyif(:Aphy == 3e4, 2150),
        @onlyif(:Aphy == 8e4, 1550),
        @onlyif(:Aphy == 3e5, 1000),
    ],
)
dicts = dict_list(parameters)

@everywhere function run_abm(config::Dict)
    @unpack mot, dt, U, λ, L, PER, Cb, α, Aphy = config
    γ = 150
    model = setup_abm_community(; mot, dt, U, λ, L, α, Aphy, PER, Cb, γ, n=100_000)
    simtime = 60 # minutes
    nsteps = round(Int, simtime * 60 / dt)
    c(a) = concentration(model)(a, model) - Cb
    adata = [x, y, z, c]
    # do not collect during first 15 minutes equilibration
    # then collect every 45 seconds
    eqtime = round(Int, 15*60 / dt)
    sampletime = round(Int, 45 / dt)
    when(model, t) = t >= eqtime && t % sampletime == 0
    adf, = run!(model, nsteps; adata, when)
    # rename because local functions don't respect names
    colnames = names(adf)
    rename!(adf,
        Symbol(colnames[end]) => :c,
    )
    adf
end

pmap(dicts) do config
    @show config
    data = run_abm(config)
    on_cluster = haskey(ENV, "SCRATCH")
    fileout = datadir(
        "comm_7classes",
        savename("comm", config, "csv");
        on_cluster
    )
    wsave(fileout, data)
end

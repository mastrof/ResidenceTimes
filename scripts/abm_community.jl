using Distributed
@everywhere using DrWatson
@everywhere @quickactivate "ResidenceTimes"
@everywhere begin
    using ResidenceTimes
    using DataFrames
    using Agents
    using MicrobeAgents
end

parameters = Dict(
    :dt => 0.1, # s
    :mot => ["RRF", "RR", "RT"],
    :U => [ # μm/s
        @onlyif(:mot == "RRF", [10, 20, 40])...,
        @onlyif(:mot == "RR", [10, 40])...,
        @onlyif(:mot == "RT", [10, 40])...,
    ],
    :λ => 2.2, # 1/s
    :PER => 0.5,
    :Cb => 0.0, # μM
    :α => 1.0,
    :Aphy => Int.([3e4, 3e5]), # cells/mL
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
    simtime = 60 # minutes
    # exposure is only accumulated after the first 15 minutes (equilibration)
    eqtime = 15 # minutes
    model = setup_abm_community(; mot, dt, U, λ, L, α, Aphy, PER, Cb, γ, n=100_000,
        teq=eqtime*60
    )
    nsteps = round(Int, simtime * 60 / dt)
    neq = model.exposure_from
    # full model state at start, end of equilibration and end of the run;
    # the fine pass replays selected bacteria from the equilibration checkpoint
    on_cluster = haskey(ENV, "SCRATCH")
    statefile(tag) = datadir(
        "comm_7classes_states",
        savename("state-$tag", config, "jls");
        on_cluster
    )
    save_checkpoint(statefile("init"), model)
    step!(model, neq)
    save_checkpoint(statefile("eq"), model)
    step!(model, nsteps - neq)
    save_checkpoint(statefile("end"), model)
    DataFrame(id=collect(allids(model)), exposure=model.exposure)
end

pmap(dicts) do config
    @show config
    data = run_abm(config)
    on_cluster = haskey(ENV, "SCRATCH")
    fileout = datadir(
        "comm_7classes_exposure",
        savename("comm", config, "csv");
        on_cluster
    )
    wsave(fileout, data)
end

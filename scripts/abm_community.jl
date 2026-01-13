using Distributed
@everywhere using DrWatson
@everywhere @quickactivate "ResidenceTimes"
@everywhere begin
    using ResidenceTimes
    using DataFrames
    using MicrobeAgents
end

mot = ["RRF"]
dt = [0.1] # s
U = [10, 40] # μm/s
λ = [2.2] # 1/s
L = [1000] # μm
PER = [0.5]
Cb = [0.03] # μM
α = [0.85]
Aphy = [3e5, 1e6, 5e6] # cells/mL
allparams = @strdict mot dt U λ L PER Cb α Aphy
dicts = dict_list(allparams)

@everywhere function run_abm(config::Dict)
    @unpack mot, dt, U, λ, L, PER, Cb, α, Aphy = config
    γ = 150
    model = setup_abm_community(; mot, dt, U, λ, L, α, Aphy, PER, Cb, γ, n=50)
    simtime = 75 # minutes
    nsteps = round(Int, simtime * 60 / dt)
    c(a) = concentration(model)(position(a), model)
    adata = [c]
    # do not collect during first 15 minutes equilibration
    # then collect every 15 seconds
    when(model, t) = t >= 0*6_000 && t % 150 == 0
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
    fileout = datadir("comm", savename("comm", config, "csv"); on_cluster)
    wsave(fileout, data)
end

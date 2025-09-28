using Distributed
@everywhere using DrWatson
@everywhere @quickactivate "ResidenceTimes"
@everywhere begin
    using ResidenceTimes
    using DataFrames
    using MicrobeAgents
end

R = [1, 3, 10] # μm
mot = ["RT", "RR", "RRF"]
dt = [0.1] # s
U = [25, 45] # μm/s
L = [1000] # μm
Cs = [0, 1, 5] # μM
allparams = @strdict R mot dt U L Cs
dicts = dict_list(allparams)

@everywhere function run_abm(config::Dict)
    @unpack R, mot, dt, U, L, Cs = config
    γ = 150
    model = setup_abm(; R, mot, dt, U, L, Cs, γ, n=50_000)
    simtime = 60 # minutes
    nsteps = round(Int, simtime * 60 / dt)
    # radial distance from source
    r(a) = distance(a, chemoattractant(model).origin, model)
    # c(a) = concentration(model)(position(a), model)
    # relative concentration value at position
    function _c(a)
        r = distance(a, chemoattractant(model).origin, model)
        exp(-(r-R)/γ) * R / r
    end
    c(a) = _c(a) # HACK: somehow necessary for correct naming
    adata = [r, c]
    # do not collect during first 30 minutes equilibration
    # then collect every 20 seconds
    when(model, t) = t >= 18_000 && t % 200 == 0
    adf, = run!(model, nsteps; adata, when) # = every 10s
    # rename because local functions don't respect names
    colnames = names(adf)
    rename!(adf,
        Symbol(colnames[end-1]) => :r,
        Symbol(colnames[end]) => :c,
    )
    adf
end

pmap(dicts) do config
    @show config
    data = run_abm(config)
    on_cluster = haskey(ENV, "SCRATCH")
    fileout = datadir("brumley", savename("abm", config, "csv"); on_cluster)
    wsave(fileout, data)
end

using Distributed
@everywhere using DrWatson
@everywhere @quickactivate "ResidenceTimes"
@everywhere begin
    using ResidenceTimes
    using DataFrames
    using MicrobeAgents
end

#== exposure distributions ==#
R = [0.5, 1, 3, 10] # μm
mot = ["RT", "RR", "RRF"]
dt = [0.1] # s
U = [10, 20, 40] # μm/s
# λ = [2.2, 0.1] # 1/s
λ = [2.2] # 1/s
L = [1000] # μm
# Cs = [0, 0.1, 1, 3] # μM
# Cs = [0.0, 1.0] # μM
Cs = [
    @onlyif("R" == 0.5, 0.24),
    @onlyif("R" == 1.0, 0.58),
    @onlyif("R" == 3.0, 2.4),
    @onlyif("R" == 10.0, 11.0),
]
Cb = [0.03] # μM
#== relative exposure vs R --- only RRF ==#
# R = [0.5] # μm
# R = [
#     2, 4, 5, 7,
#     @onlyif("U" == 60 || "U" == 80, [1, 3, 10])...
# ] # μm
# mot = ["RRF"]
# dt = [0.1] # s
# U = [10, 20, 40, 60, 80] # μm/s
# λ = [2.2, 0.1] # 1/s
# L = [1000] # μm
# Cs = [0, 1] # μM
# Cb = [0.03] # μM
#== diffusivity matching RR vs RT ==#
# R = [1.0, 3.0] # μm
# mot = ["RR"]
# dt = [0.1] # s
# λ = [2.2] # 1/s
# L = [1000] # μm
# Cs = [1] # μM
# Cb = [0.03] # μM
# Dr = 0.035 # rad²/s --- same in all simulations
# τ = 1/λ[1] # s
# U = [10, 20, 40] ./ sqrt((1+2*Dr*τ)/(2*(1+Dr*τ))) # μm/s
allparams = @strdict R mot dt U λ L Cs Cb
dicts = dict_list(allparams)

@everywhere function run_abm(config::Dict)
    @unpack R, mot, dt, U, λ, L, Cs, Cb = config
    γ = 150
    model = setup_abm(; R, mot, dt, U, λ, L, Cs, Cb, γ, n=50_000)
    simtime = 75 # minutes
    nsteps = round(Int, simtime * 60 / dt)
    # radial distance from source
    r(a) = distance(a, chemoattractant(model).origin, model)
    # normalized concentration value at position
    cellradius = radius(model[1])
    c_max = exp(-cellradius/γ) * R / (R+cellradius)
    function _c(a)
        r = distance(a, chemoattractant(model).origin, model)
        exp(-(r-R)/γ) * R / r / c_max
    end # =1 when cell and phytoplankton are in contact
    c(a) = _c(a) # HACK: somehow necessary for correct naming
    adata = [r, c]
    # do not collect during first 30 minutes equilibration
    # then collect every 45 seconds
    when(model, t) = t >= 18_000 && t % 450 == 0
    adf, = run!(model, nsteps; adata, when)
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
    fileout = datadir("abm", savename("abm", config, "csv"); on_cluster)
    wsave(fileout, data)
end

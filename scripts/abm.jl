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
Cs = [0, 1, 10] # μM
allparams = @strdict R mot dt U L Cs
dicts = dict_list(allparams)

@everywhere function run_abm(config::Dict)
    @unpack R, mot, dt, U, L, Cs = config
    γ = 150
    model = setup_abm(; R, mot, dt, U, L, Cs, γ, n=50_000)
    simtime = 20 # minutes
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
    adf, = run!(model, nsteps; adata, when=100) # = every 10s
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
    fileout = datadir("brumley", savename(config, "csv"); on_cluster)
    wsave(fileout, data)
end


# let Δ=5
#     fig = Figure(size=(800,800))
#     ax = Axis(fig[1,1])
#     ylims!(ax, 0.5, 3)
#     display(fig)
#     t = Observable(0)
#     timestamp = @lift(string($t))
#     text!(ax, 0.02, 0.95; text=timestamp, space=:relative)
#     r = range(R, 600; length=120)
#     bandwidth = 20
#     k0 = pdf(kde(adf[adf.time .== 0, :radial_distance]; bandwidth), r)
#     kt = @lift(
#         pdf(kde(
#             adf[$t-Δ*when .<= adf.time .<= $t, :radial_distance];
#             bandwidth
#         ), r)
#     )
#     g = @lift($kt ./ k0)
#     lines!(ax, r, g; linewidth=3)
#     for s in unique(adf.time)[1:Δ:end]
#         sleep(1/10)
#         t[] = s
#     end
# end

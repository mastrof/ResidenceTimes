using DrWatson
@quickactivate "ResidenceTimes"
push!(LOAD_PATH, projectdir("analysis")) # Makie lives in the analysis workspace project
using ResidenceTimes
using CSV
using DataFrames
using LinearAlgebra
using StatsBase
using KernelDensity
using LaTeXStrings
using GLMakie
set_theme!(theme_latexfonts();
    fontsize=18,
    Axis=(
        xgridvisible=false,
        ygridvisible=false,
        xticksize=-5,
        yticksize=-5,
        xticksmirrored=true,
        yticksmirrored=true,
        titlesize=18,
    ),
    Legend=(
        framecolor=:transparent,
        backgroundcolor=RGBAf(1,1,1,0.6),
        patchsize=(40, 16),
        titlesize=18,
        titlefont=:regular,
        labelsize=18,
    ),
    Lines=(
        linewidth=6,
    )
)
it(s) = rich(s; font=:italic)

unique_phyto(df::AbstractDataFrame) = unique_phyto(groupby(df, :id))
unique_phyto(gdf::GroupedDataFrame) = [unique_phyto(g) for g in gdf]
function unique_phyto(g::SubDataFrame)
    length(unique(g.idp))
end

filenames = filter(
    s -> startswith(s, "dist"),
    readdir(datadir("comm_7classes_proc"))
)[1:2]
comm = let
    ks = ["Aphy", "L", "α"]
    parameters = parse_savename(filenames[1])[2]
    comm_parameters = Dict(ks .=> get.(Ref(parameters), ks, nothing))
    CSV.read(datadir(savename("phytoplankton", comm_parameters, "csv")), DataFrame)
end

data = map(filenames) do fname
    df = CSV.read(
        datadir("comm_7classes_proc", fname),
        DataFrame;
    )
    ids, es = eachcol(CSV.read(
        datadir(
            "comm_7classes_ids",
            replace(fname, "dist"=>"ids1-99")
        ),
        DataFrame
    ))
    q = sortperm(es)
    lo = q[1:1000]
    hi = q[1001:end]
    Dict(
        :lo => subset(df, :id => ByRow(id -> id ∈ ids[lo])),
        :hi => subset(df, :id => ByRow(id -> id ∈ ids[hi])),
    )
end
df_lowU = data[1]
df_hiU = data[2]

let
    colors = [
        "#B133DE",
        "#205C80",
    ]
    fig = Figure(size=(300,500))
    Rs = unique(comm.R)
    let
        x = -0.5:21.5
        ax = Axis(fig[1,1];
            xlabel=rich("Unique phytoplankton"),
            ylabel="probability density",
            yticklabelsvisible=false,
        )
        y = kde(unique_phyto(df_lowU[:hi]))
        lines!(ax, x, x->pdf(y,x);
            color=colors[1],
            linewidth=4,
            label="best",
        )
        y = kde(unique_phyto(df_lowU[:lo]))
        lines!(ax, x, x->pdf(y,x);
            color=colors[1],
            linewidth=4,
            linestyle=:dash,
            label="worst",
        )
        axislegend(ax, "10 μm/s"; position=:rt, patchsize=(40,16))
    end
    let
        x = -2:4:102
        ax = Axis(fig[2,1];
            xlabel=rich("Unique phytoplankton"),
            ylabel="probability density",
            yticklabelsvisible=false,
        )
        y = kde(unique_phyto(df_hiU[:hi]))
        lines!(ax, x, x->pdf(y,x);
            color=colors[2],
            linewidth=4,
            label="best"
        )
        y = kde(unique_phyto(df_hiU[:lo]))
        lines!(ax, x, x->pdf(y,x);
            color=colors[2],
            linewidth=4,
            linestyle=:dash,
            label="worst",
        )
        axislegend(ax, "40 μm/s"; position=:lt, patchsize=(40,16))
    end
    fig
end

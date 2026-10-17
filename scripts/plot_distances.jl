using DrWatson
@quickactivate "ResidenceTimes"
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

unique_phyto(gdf::GroupedDataFrame) = [unique_phyto(g) for g in gdf]
function unique_phyto(g::SubDataFrame)
    length(unique(g.idp))
end

function meandist(df::AbstractDataFrame)
    gdf = groupby(df, :id)
    [mean(g.d) for g in gdf]
end

filenames = filter(
    s -> startswith(s, "dist"),
    readdir(datadir("comm_7classes_proc"))
)

data = map(filenames[1:2]) do fname
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
    fig = Figure(size=(600,500))
    x = range(0, 420; length=2500)
    let
        ax = Axis(fig[1:2,2:3];
            xlabel=rich(it("d"), " (μm)"),
            ylabel="probability density",
            yticklabelsvisible=false,
        )
        ax_b1 = Axis(fig[1,1];
            ylabel=rich(it("d"), " (μm)"),
        )
        ax_b2 = Axis(fig[2,1];
            xlabel=rich(it("t"), " (min)"),
            ylabel=rich(it("d"), " (μm)"),
        )
        y = kde(meandist(df_lowU[:hi]))
        lines!(ax, x, x -> pdf(y,x);
            color=colors[1],
            label="best"
        )
        gdf = groupby(df_lowU[:hi], :id)
        for i in [2,3,5]
            t = (gdf[i].time .- gdf[i].time[1]) .* (0.1/60)
            d = gdf[i].d
            s = 6 # subsample
            lines!(ax_b1, t[1:s:end], d[1:s:end]; linewidth=1)
        end
        y = kde(meandist(df_lowU[:lo]))
        lines!(ax, x, x -> pdf(y,x);
            color=colors[1],
            label="worst",
            linewidth=4,
            linestyle=:dash
        )
        gdf = groupby(df_lowU[:lo], :id)
        for i in [1, 2, 3]
            t = (gdf[i].time .- gdf[i].time[1]) .* (0.1/60)
            d = gdf[i].d
            s = 6 # subsample
            lines!(ax_b2, t[1:s:end], d[1:s:end]; linewidth=1)
        end
        axislegend(ax, "10 μm/s"; position=:rt, patchsize=(40,16))
    end
    let
        ax = Axis(fig[3:4,2:3];
            xlabel=rich(it("d"), " (μm)"),
            ylabel="probability density",
            yticklabelsvisible=false,
        )
        ax_b1 = Axis(fig[3,1];
            ylabel=rich(it("d"), " (μm)"),
        )
        ax_b2 = Axis(fig[4,1];
            xlabel=rich(it("t"), " (min)"),
            ylabel=rich(it("d"), " (μm)"),
        )
        y = kde(meandist(df_hiU[:hi]))
        lines!(ax, x, x -> pdf(y,x);
            color=colors[2],
            label="best"
        )
        gdf = groupby(df_hiU[:hi], :id)
        for i in [8,5,7]
            t = (gdf[i].time .- gdf[i].time[1]) .* (0.1/60)
            d = gdf[i].d
            s = 6 # subsample
            lines!(ax_b1, t[1:s:end], d[1:s:end]; linewidth=1)
        end
        y = kde(meandist(df_hiU[:lo]))
        lines!(ax, x, x -> pdf(y,x);
            color=colors[2],
            label="worst",
            linewidth=4,
            linestyle=:dash
        )
        gdf = groupby(df_hiU[:lo], :id)
        # for i in [8,3,5] # oligo
        for i in [10,1,9] # eu
            t = (gdf[i].time .- gdf[i].time[1]) .* (0.1/60)
            d = gdf[i].d
            s = 6 # subsample
            lines!(ax_b2, t[1:s:end], d[1:s:end]; linewidth=1)
        end
        axislegend(ax, "40 μm/s"; position=:rt, patchsize=(40,16))
    end
    fig
end

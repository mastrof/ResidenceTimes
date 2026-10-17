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
    fontsize=22,
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

filenames = filter(
    s -> startswith(s, "gyration"),
    readdir(datadir("comm_7classes_proc"))
)

data = map(filenames[3:4]) do fname
    df = CSV.read(
        datadir("comm_7classes_proc", fname),
        DataFrame;
    )
    ids, es = eachcol(CSV.read(
        datadir(
            "comm_7classes_ids",
            replace(fname, "gyration"=>"ids1-99")
        ),
        DataFrame
    ))
    q = sortperm(es)
    lo = q[1:1000]
    hi = q[1001:end]
    Dict(:lo => df[lo,:], :hi => df[hi,:])
end

df_lowU = data[1]
df_hiU = data[2]

let
    colors = [
        "#B133DE",
        "#205C80",
    ]
    fig = Figure(size=(467,500))
    # x = range(0.05, 1; length=1000)
    # bandwidth=4e-2
    # Rg0(U, N=27000, λ=2.2) = sqrt(N/6) * (U/λ)
    let
        ax = Axis(fig[1,1];
            # xlabel="Gyration radius (μm)",
            ylabel="probability density",
        )
        x = range(20, 200; length=1000)
        # y = kde(log10.(df_lowU[:hi].V); bandwidth)
        y = kde(df_lowU[:hi].Rg)
        lines!(ax, x, x -> pdf(y,x);
            color=colors[1],
            label="best"
        )
        vlines!(ax, mean(df_lowU[:hi].Rg); color=colors[1])
        # y = kde(log10.(df_lowU[:lo].V); bandwidth)
        y = kde(df_lowU[:lo].Rg)
        lines!(ax, x, x -> pdf(y,x);
            color=colors[1],
            label="worst",
            linewidth=4,
            linestyle=:dash
        )
        vlines!(ax, mean(df_lowU[:lo].Rg); color=colors[1], linestyle=:dash)
        axislegend(ax, "10 μm/s"; position=:rt, patchsize=(40,16))
        text!(ax, 0.02, 0.8; text="A",
            space=:relative,
            font=:bold,
            fontsize=28
        )
        # vlines!(ax, Rg0(10); color=:black)
        println(mean(df_lowU[:hi].Rg) - mean(df_lowU[:lo].Rg))
    end
    let
        ax = Axis(fig[2,1];
            xlabel=rich("Gyration radius, ", it("R"), subscript("g"), " (μm)"),
            # xlabel=LaTeXString("{\\mathrm{Rescaled gyration radius, }} \\textit{R_g/R_g^0}"),
            ylabel="probability density",
        )
        x = range(150, 850; length=1000)
        # y = kde(log10.(df_hiU[:hi].V); bandwidth)
        y = kde(df_hiU[:hi].Rg)
        lines!(ax, x, x -> pdf(y,x);
            color=colors[2],
            label="best"
        )
        vlines!(ax, mean(df_hiU[:hi].Rg); color=colors[2])
        # y = kde(log10.(df_hiU[:lo].V); bandwidth)
        y = kde(df_hiU[:lo].Rg)
        lines!(ax, x, x -> pdf(y,x);
            color=colors[2],
            label="worst",
            linewidth=4,
            linestyle=:dash
        )
        vlines!(ax, mean(df_hiU[:lo].Rg); color=colors[2], linestyle=:dash)
        axislegend(ax, "40 μm/s"; position=:rt, patchsize=(40,16))
        text!(ax, 0.02, 0.8; text="B",
            space=:relative,
            font=:bold,
            fontsize=28
        )
        println(mean(df_hiU[:hi].Rg) - mean(df_hiU[:lo].Rg))
    end
    fig
end

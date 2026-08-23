using DrWatson
@quickactivate "ResidenceTimes"
using ResidenceTimes
using CSV
using DataFrames
using StatsBase
using Distributions
using KernelDensity
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
        titlefont=:italic,
        labelsize=18,
    ),
    Lines=(
        linewidth=6,
    )
)

filenames = filter(
    s -> startswith(s, "encounters"),
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
            replace(fname, "encounters"=>"ids1-99")
        ),
        DataFrame
    ))
    q = sortperm(es)
    lo = q[1:1000]
    hi = q[1001:end]
    (
        Dict(
            :lo => df.encounters[lo],
            :hi => df.encounters[hi],
        ),
        Dict(
            :lo => df.residencetime[lo],
            :hi => df.residencetime[hi],
        ),
    )
end
enc_lowU, resperenc_lowU = data[1]
enc_hiU, resperenc_hiU = data[2]

# encounters low vs high speed
let
    colors = [
        "#B133DE",
        "#205C80",
    ]
    fig = Figure(size=(467,500))
    let
        ax = Axis(fig[1,1];
            xlabel="hotspot encounter rate (1/hr)",
            ylabel="probability density",
        )
        y1 = enc_lowU[:hi] .* (60/45)
        y2 = enc_hiU[:hi] .* (60/45)
        hist!(ax, y1;
            bins=0:10:400,
            normalization=:pdf,
            color=colors[1]*"77",
        )
        hist!(ax, y2;
            bins=0:10:400,
            normalization=:pdf,
            color=colors[2]*"77",
        )
        f1 = fit(Gamma, y1)
        f2 = fit(Gamma, y2)
        x = 0..400
        lines!(ax, x, x -> pdf(f1, x);
            color=colors[1],
            # linewidth=3,
            label="10 μm/s",
        )
        lines!(ax, x, x -> pdf(f2, x);
            color=colors[2],
            # linewidth=3,
            label="40 μm/s",
        )
        axislegend(ax; position=:ct, patchsize=(40,16))
        ylims!(ax, 0, nothing)
        xlims!(ax, 0, nothing)
        text!(ax, 0.9, 0.8; text="A",
            space=:relative,
            font=:bold,
            fontsize=28
        )
    end
    let
        ax = Axis(fig[2,1];
            xlabel="mean encounter duration (s)",
            ylabel="probability density",
            xscale=log10,
            yscale=log10,
        )
        y1 = resperenc_lowU[:hi] .* 0.1
        y2 = resperenc_hiU[:hi] .* 0.1
        hist!(ax, y1;
            bins=logrange(1, 1e3; length=30),
            normalization=:pdf,
            color=colors[1]*"77",
        )
        hist!(ax, y2;
            bins=logrange(1, 1e3; length=30),
            normalization=:pdf,
            color=colors[2]*"77",
        )
        ylims!(ax, 1e-4, nothing)
        xlims!(ax, 1, nothing)
        text!(ax, 0.9, 0.8; text="B",
            space=:relative,
            font=:bold,
            fontsize=28
        )
    end
    fig
end

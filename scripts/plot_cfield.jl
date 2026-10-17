using DrWatson
@quickactivate "ResidenceTimes"
using ResidenceTimes
using ResidenceTimes: leaked_concentration, leakage_rate
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
        linewidth=4,
    )
)
it(s) = rich(s; font=:italic)

dist(a, b) = hypot((a.-b)...)

conc = ResidenceTimes._concentration
γ = 150.0

filenames = filter(startswith("phytoplankton"), readdir(datadir()))
fname = filenames[1]
comm = CSV.read(datadir(fname), DataFrame)[1:301,:]
transform!(comm,
    :R => ByRow(R -> leaked_concentration(R, 0.3)) => :Cs,
    :R => ByRow(R -> leakage_rate(R, 0.3)) => :L,
)

let
    fig = Figure(size=(700,400))
    ax1 = Axis(fig[1:3,1:3];
        xlabel=rich(it("r"), " (μm)"),
        ylabel=rich(it("C(r)"), " (μM)"),
    )
    r = range(0, 75; length=1000)
    Rs = unique(comm.R)
    for R in Rs
        Cs = leaked_concentration(R, 0.3)
        c = conc.(r[r.>R], R, Cs, γ)
        lines!(ax1, r[r.>R], c; label=string(R))
    end
    axislegend(ax1, "radius (μm)"; position=:rt)
    #
    ax2 = Axis(fig[1,4:6];
        xlabel=rich("radius (μm)"),
        ylabel=rich("cells/mL"),
        xticks=1:2:7,
        yticks=3:5,
        ytickformat=(ys -> map(y -> rich("10", superscript(string(Int(y)))), ys))
    )
    V = parse_savename(fname)[2]["L"]^3 # μm^3
    V *= 1e-12 # μm^3 → mL
    N = [count(==(R), comm.R) for R in Rs] ./ V
    barplot!(ax2, Rs, log10.(N);
        color=:grey,
        width=Rs,
        fillto=2.8,
    )
    ylims!(ax2, 2.8, nothing)
    #
    ax3 = Axis(fig[2:3,4:6])
    hidedecorations!(ax3)
    hidespines!(ax3)
    xs = (0, 800) .+ 20
    ys = (0, 600) .+ 0
    zs = (0, 300)
    subcomm = subset(comm,
        :x => ByRow(x -> xs[1]<x<xs[2]),
        :y => ByRow(y -> ys[1]<y<ys[2]),
        :z => ByRow(z -> zs[1]<z<zs[2]),
    )
    δ = 5
    xrng = range(xs...; step=δ)
    yrng = range(ys...; step=δ)
    z = (zs[1]+zs[2])/2
    f = zeros(length(xrng), length(yrng))
    for CI in CartesianIndices(f)
        i, j = CI.I
        x = xrng[i]
        y = yrng[j]
        for phyto in eachrow(subcomm)
            d = dist([x,y,z], [phyto.x, phyto.y, phyto.z])
            R = phyto.R
            Cs = phyto.Cs
            f[CI] += conc(d, R, Cs, γ)
        end
    end
    hm = heatmap!(ax3, xrng, yrng, f.*1e3)
    Colorbar(fig[2:3,7], hm;
        label="nM"
    )
    fig
end

let
    mots = ["RRF", "RR", "RT"]
    fig = Figure(size=(360,180))
    ax = Axis(fig[1,1];
        ylabel=rich(it("D"), " (μm²/s)"),
        xticks=1:3,
        xtickformat=(ys -> [mots[Int(y)] for y in ys])
    )
    τrun = (1 / 2.2)u"s"
    τrot = (1 / (2*0.035))u"s"
    τeff = τrun*τrot/(τrun+τrot)
    v = [10, 20, 40]u"μm/s"
    D_RT = @. v^2 * τeff / 3 |> u"μm^2/s" |> ustrip
    D_RR = @. v^2 * τeff / 6 |> u"μm^2/s" |> ustrip
    D_RRF = @. v^2/6 * (1/τrun + 2/τrot) / (1/τrun + 1/τrot)^2 |> u"μm^2/s" |> ustrip
    markersize = 16
    marker = [:circle, :rect, :diamond]
    scatter!(ax, fill(1, length(v)), D_RRF; markersize, marker)
    scatter!(ax, fill(2, length(v)), D_RR; markersize, marker)
    scatter!(ax, fill(3, length(v)), D_RT; markersize, marker)
    Legend(fig[1,1], [
            MarkerElement(; color=:black, marker=marker[1], markersize),
            MarkerElement(; color=:black, marker=marker[2], markersize),
            MarkerElement(; color=:black, marker=marker[3], markersize),
        ],
        ["10 μm/s", "20 μm/s", "40 μm/s"];
        valign=:top,
        halign=:left,
        tellwidth=false,
        tellheight=false,
        backgroundcolor=:transparent,
        patchlabelgap=0,
    )
    fig
end

using DrWatson
@quickactivate "ResidenceTimes"
push!(LOAD_PATH, projectdir("analysis")) # Makie lives in the analysis workspace project
using ResidenceTimes
conc = ResidenceTimes._concentration
const γ = 150.0
const Rmin = 0.65
const Rmax = 10.0
using IterTools: partition
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

filenames = filter(startswith("phytoplankton"), readdir(datadir()))
fname = filenames[1]
comm = CSV.read(datadir(fname), DataFrame)#[1:301,:]
field = load_community_field(datadir(replace(fname, "phytoplankton"=>"cfield", "csv"=>"jld2")))

let
    fig = Figure(size=(700,400))
    ax1 = Axis(fig[1:3,1:3];
        xlabel=rich(it("r"), " (μm)"),
        ylabel=rich(it("C(r)"), " (nM)"),
    )
    r = range(0, 75; length=1000)
    Rs = unique(comm.R)
    Cs = unique(comm.Cs)
    for i in eachindex(Rs)
        R = Rs[i]
        C = Cs[i]
        c = conc.(r[r.>R], R, C, γ) .* 1e3
        lines!(ax1, r[r.>R], c; label=string(R))
    end
    axislegend(ax1, "Cell radius (μm)"; position=:rt)
    #
    ax2 = Axis(fig[1,4:6];
        xlabel=rich("Cell radius (μm)"),
        ylabel=rich("cells/mL"),
        xticks=1:2:7,
        yticks=3:5,
        ytickformat=(ys -> map(y -> rich("10", superscript(string(Int(y)))), ys))
    )
    V = parse_savename(fname)[2]["L"]^3 # μm^3
    V *= 1e-12 # μm^3 → mL
    N = [count(==(R), comm[.~comm.image, :R]) for R in Rs] ./ V
    Rbins = reverse(collect(partition(logrange(Rmin, Rmax; length=8), 2, 1)))
    for i in eachindex(N)
        r1, r2 = Rbins[i+1]
        dr = r2-r1
        r = (r1+r2)/2
        n = N[i]
        barplot!(ax2, [r], [log10(n)];
            color=Makie.wong_colors()[i],
            width=2dr,
            fillto=2.8
        )
    end
    ylims!(ax2, 2.8, nothing)
    #
    ax3 = Axis(fig[2:3,4:6])
    hidedecorations!(ax3)
    hidespines!(ax3)
    xs = (0, 960) .+ 20
    ys = (0, 720) .+ 0
    zs = (0, 300)
    δ = 2
    xrng = range(xs...; step=δ)
    yrng = range(ys...; step=δ)
    z = (zs[1]+zs[2])/2
    f = zeros(length(xrng), length(yrng))
    for CI in CartesianIndices(f)
        i, j = CI.I
        x = xrng[i]
        y = yrng[j]
        f[CI] = field(SVector(x, y, z))
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

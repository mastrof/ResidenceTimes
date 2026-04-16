using DrWatson
using ResidenceTimes
using CSV
using DataFrames
using KernelDensity
using LaTeXStrings
using Colors
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
        labelsize=18
    ),
    Lines=(
        linewidth=4,
    )
)

# exposure normalization
function E0(R; ntimepoints=62, γ=150, L=1e3)
    # a = 4π/3 * R^3 # volume of the cell
    # b = (R - γ)*(γ/L * exp(-(L-R)/γ) - γ/R)
    # (a + b) * ntimepoints
    rmin = 0.5
    c = (R - γ) * (γ/L * exp(-(L-R)/γ) - γ/rmin * exp(-(rmin-R)/γ))
end

function plot_exposurehist!(ax, exposure, R, Cs, mot, λs, Us)
    colors = cgrad(:Paired_10)
    for i in eachindex(exposure)
        label = join(string.([Us[i], λs[i]]), ", ")
        color = colors[i]
        x = range(-3.5, 2.5; length=400)
        e = log10.(exposure[i])
        # x = range(-8.5, -1.5; length=400)
        # e = log10.(exposure[i] ./ E0(R))
        linestyle = iseven(i) ? :solid : Linestyle([0,4,5])
        linewidth = iseven(i) ? 4 : 3
        lines!(ax, x, pdf(kde(e; bandwidth=0.09), x);
            label, color, linestyle, linewidth
        )
    end
end

fig = Figure(size=(1200,800))

#==
A row with 3 panels.
Each panel shows the exposure pdf for RR motility at a given R.
Each one includes different U and λ values at fixed Cs=1.0.
==#
g_top = fig[1,:] = GridLayout()
Rs = [1, 2, 3]
ax_top = [Axis(g_top[1,i]) for i in 1:3]
for (i,R) in enumerate(Rs)
    Cs = 1.0
    mot = "RRF"
    datasets = collect_datasets(datadir("exposure");
        R=R,
        Cs=Cs,
        mot=mot,
    )
    filter!(:λ => l -> l .!= 0.5, datasets)
    filter!(:U => u -> u .<= 40, datasets) # [10, 20, 40]
    sort!(datasets, [:U, :λ])
    gdata = groupby(datasets, [:U, :λ])
    exposure = [g.c_sum for g in gdata]
    Us = [first(g.U) for g in gdata]
    λs = [first(g.λ) for g in gdata]
    plot_exposurehist!(ax_top[i], exposure, R, Cs, mot, λs, Us)
    ylims!(ax_top[i], 1e-4, 2)
    ax_top[i].yscale = log10
    ax_top[i].title = LaTeXString(join(
        [
            LaTeXString("\$R\$=$(R)\\,μm"),
            LaTeXString("$(mot) motility"),
            LaTeXString("\$C_\\mathrm{S}\$=$(Cs)\\,μM"),
        ],
        LaTeXString(",\\;")
    ))
    ax_top[i].xlabel = "log exposure"
end
ax_top[2].yticklabelsvisible = false
ax_top[3].yticklabelsvisible = false
axislegend(ax_top[3],
    LaTeXString("\$U\\,\\mathrm{(μm/s)},\\;λ\\,\\mathrm{(1/s)}\$");
    position=:rt
)
ax_top[1].ylabel = "PDF"

#==
Row with 2 panels
1. Average exposure vs R for RRF motility at Cs=1
2. Exposure distributions for RR/RRF/RRT at R=1, Cs=1, U=20, λ=2.2
==#
g_mid = fig[2,:] = GridLayout()
ax_mid = [Axis(g_mid[1,i]) for i in 1:2]
mots = ["RR", "RRF", "RT"]
Rs = [1, 2, 3, 4, 5, 7, 10]
let mot="RRF", Cs=2.0
    datasets = collect_datasets(datadir("exposure");
        Cs=Cs,
        mot=mot,
    )
    filter!(:λ => l -> l .!= 0.5, datasets) # [0.1, 2.2]
    sort!(datasets, [:R, :U, :λ])
    gdata = groupby(datasets, [:U, :λ])
    exposure = map(collect(gdata)) do dataset
        gdfR = groupby(dataset, :R)
        [mean(gR.c_sum) for gR in gdfR]# ./ E0.(Rs)
    end
    Us = [first(g.U) for g in gdata]
    λs = [first(g.λ) for g in gdata]
    datasets0 = collect_datasets(datadir("exposure");
        Cs=0,
        mot=mot,
    )
    sort!(datasets0, [:R, :U, :λ])
    gdata0 = groupby(datasets0, [:U, :λ])
    exposure0 = map(collect(gdata0)) do dataset0
        gdfR = groupby(dataset0, :R)
        [mean(gR.c_sum) for gR in gdfR]# ./ E0.(Rs)
    end
    colors = cgrad(:Paired_10)
    markers = [
        fill(:circle, 2);
        fill(:diamond, 2);
        fill(:utriangle, 2);
        fill(:rect, 2);
    ]
    for j in eachindex(exposure)
        label = join(string.([Us[j], λs[j]]), ", ")
        color = colors[j]
        marker = markers[j]
        scatterlines!(ax_mid[1], Rs, exposure[j] ./ (0 .+ 1exposure0[j]);
            label, color, linewidth=4, marker, markersize=18
        )
    end
    ax_mid[1].title = LaTeXString(join(
        [
            LaTeXString("$(mot) motility"),
            LaTeXString("\$C_\\mathrm{S}\$=$(Cs)\\,μM")
        ],
        LaTeXString(",\\;")
    ))
end
axislegend(ax_mid[1],
    LaTeXString("\$U\\,\\mathrm{(μm/s)},\\;λ\\,\\mathrm{(1/s)}\$");
    position=:rt, nbanks=2
)
ax_mid[1].yscale = log10
ax_mid[1].xlabel = LaTeXString("\$R\\,\\mathrm{(μm)}\$")
ax_mid[1].ylabel = LaTeXString("Average exposure gain")
ax_mid[1].yticks = [1, 3, 10, 30, 100]

let Cs=1.0, λ=2.2, U=20, R=1
    for (i,mot) in enumerate(mots)
        datasets = collect_datasets(datadir("exposure");
            Cs=Cs,
            R=R,
            λ=λ,
            U=U,
            mot=mot,
        )
        gdata = groupby(datasets, [:U])
        exposure = [g.c_sum for g in gdata]
        colors = [
            cgrad(:Blues)[3:2:7];
            cgrad(:Greens)[3:2:7];
            cgrad(:Reds)[3:2:7]
        ]
        for j in eachindex(exposure)
            label = LaTeXString("$(mot)")
            color = colors[3*(i-1)+j+1]
            x = range(-3.5, 2.5; length=400)
            e = log10.(exposure[j])
            lines!(ax_mid[2], x, pdf(kde(e; bandwidth=0.09), x);
                label, color,
            )
        end
    end
    ylims!(ax_mid[2], 1e-4, 2)
    ax_mid[2].yscale = log10
    ax_mid[2].title = LaTeXString(
        "\$R=$(R)\\,\\mathrm{μm},\\,U=$(U)\\,\\mathrm{μm/s},\\,λ=$(λ)\\,\\mathrm{s^{-1}}\$"
    )
    ax_mid[2].title = LaTeXString(join(
        [
            LaTeXString("\$R\$=$(R)\\,μm"),
            LaTeXString("\$U\$=$(U)\\,μm/s"),
            LaTeXString("\$λ\$=$(λ)\\,s^{-1}"),
            LaTeXString("\$C_\\mathrm{S}\$=$(Cs)\\,μM"),
        ],
        LaTeXString(",\\;")
    ))
    ax_mid[2].xlabel = "log exposure"
    ax_mid[2].ylabel = "PDF"
end
axislegend(ax_mid[2],
    LaTeXString("\\mathrm{motility}");
    position=:rt,
)

fig

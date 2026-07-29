using DrWatson
using ResidenceTimes
using CSV
using DataFrames
using KernelDensity
using LaTeXStrings
using Printf
using StatsBase
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
    # colors = cgrad(:Paired_10)
    colors = [
        "#B133DE",
        "#D57D19",
        "#205C80",
    ]
    for i in eachindex(exposure)
        # label = join(string.([Us[i], λs[i]]), ", ")
        label = string(Us[i])
        color = colors[i]
        x = range(-2.5, 4.0; length=400)
        e = log10.(exposure[i])
        # x = range(-8.5, -1.5; length=400)
        # e = log10.(exposure[i] ./ E0(R))
        # linestyle = iseven(i) ? :solid : Linestyle([0,4,5])
        # linewidth = iseven(i) ? 4 : 3
        linestyle = :solid
        linewidth = 6
        k = kde(e; bandwidth=0.06)
        lines!(ax, x, pdf(k, x);
            label, color, linestyle, linewidth
        )
        e_avg = mean(exposure[i])
        text!(ax, 0.85, 0.35-0.075*(i-1);
            text=string(round(e_avg; sigdigits=2)),
            space=:relative,
            fontsize=18,
            color,
        )
        @printf "/** %s μm/s **/\n" string(Us[i])
        @printf "top 0.1%%: %f\n" percentile(exposure[i], 99.9)
    end
end

fig = Figure(size=(1400,400))
N = 61 # num time samples per dataset

#==
A row with 3 panels.
Each panel shows the exposure pdf for RR motility at a given R.
Each one includes different U and λ values at fixed Cs=1.0.
==#
g_top = fig[1,:] = GridLayout()
Rs = [1.0, 3, 10]
ax_top = [Axis(g_top[1,i]) for i in 1:3]
for (i,R) in enumerate(Rs)
    # Cs = 1.0
    Cs = round(ResidenceTimes.leaked_concentration(R, 0.5), sigdigits=2)
    mot = "RRF"
    datasets = collect_datasets(datadir("exposure");
        R=R,
        Cs=Cs,
        mot=mot,
        λ=2.2,
    )
    # filter!(:λ => l -> l .!= 0.5, datasets)
    filter!(:U => u -> u .<= 40, datasets) # [10, 20, 40]
    sort!(datasets, [:U, :λ])
    gdata = groupby(datasets, [:U, :λ])
    exposure = [g.c_sum ./ N .* 1e3 for g in gdata]
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
    ax_top[i].xlabel = L"\langle\,\tilde{E}\,\rangle\,\times\,1000"
    ax_top[i].xtickformat=(xs -> [LaTeXString("10^{$(Int(x))}") for x in xs])
    panel = string(('A':'C')[i])
    text!(ax_top[i], 0.03, 0.88; text=panel,
        space=:relative,
        font=:bold,
        fontsize=28
    )
end
ax_top[2].yticklabelsvisible = false
ax_top[3].yticklabelsvisible = false
axislegend(ax_top[1],
    # LaTeXString("\$U\\,\\mathrm{(μm/s)},\\;λ\\,\\mathrm{(1/s)}\$");
    LaTeXString("\$U\\,\\mathrm{(μm/s)}\$");
    position=:rt
)
ax_top[1].ylabel = "Probability density"
for ax in ax_top
    text!(ax, 0.83, 0.35+0.075;
        text="Mean:",
        space=:relative,
        fontsize=18,
        color=:black,
    )
end

fig

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

function plot_exposurehist!(ax, exposure, mot)
    # colors = cgrad(:Paired_10)
    colors = cgrad(:tab10)
    for i in eachindex(exposure)
        label = string(mot[i])
        color = colors[i]
        x = range(-2.5, 4.0; length=400)
        e = log10.(exposure[i])
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
    end
end

fig = Figure(size=(467,400))
N = 61 # num time samples per dataset

ax = Axis(fig[1,1])
let
    U = 20
    R = 3
    Cs = round(ResidenceTimes.leaked_concentration(R, 0.5), sigdigits=2)
    datasets = collect_datasets(datadir("exposure");
        R=R,
        Cs=Cs,
        U=U,
        λ=2.2,
    )
    sort!(datasets, [:U, :λ])
    gdata = groupby(datasets, [:mot])
    exposure = [g.c_sum ./ N .* 1e3 for g in gdata]
    mots = [first(g.mot) for g in gdata]
    plot_exposurehist!(ax, exposure, mots)
    ylims!(ax, 1e-4, 2)
    ax.yscale = log10
    ax.title = LaTeXString(join(
        [
            LaTeXString("\$R\$=$(R)\\,μm"),
            LaTeXString("\$U\$=$(U)\\,μm/s"),
            LaTeXString("\$C_\\mathrm{S}\$=$(Cs)\\,μM"),
        ],
        LaTeXString(",\\;")
    ))
    ax.xlabel = L"\langle\,\tilde{E}\,\rangle\,\times\,1000"
    ax.xtickformat=(xs -> [LaTeXString("10^{$(Int(x))}") for x in xs])
    # panel = string(('A':'C')[i])
    # text!(ax_top[i], 0.03, 0.88; text=panel,
    #     space=:relative,
    #     font=:bold,
    #     fontsize=28
    # )
end
axislegend(ax,
    LaTeXString("Motility");
    position=:rt
)
ax.ylabel = "Probability density"
text!(ax, 0.83, 0.35+0.075;
    text="Mean:",
    space=:relative,
    fontsize=18,
    color=:black,
)

fig

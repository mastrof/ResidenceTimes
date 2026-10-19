export setup_abm, setup_abm_community, setup_abm_community_Δche
export generate_community

function setup_abm(;
    n=10_000, # number of bacteria
    L=1e3, # domain size (μm)
    periodic=true, # periodic boundary conditions (true / false)
    R=1.0, # source radius (μm)
    Cb=0.0, # background chemoattractant (μM)
    Cs=1.0, # source chemoattractant (μM)
    γ=150.0, # chemoattractant decay lengthscale (μm)
    mot="RT", # motility pattern (RT / RR / RRF)
    dt=0.02, # integration timestep (s)
    U=45.0, # swimming speed (μm/s)
    λ=2.2, # reorientation rate (1/s)
    Drot=0.035, # rotational diffusivity (rad^2/s)
    Π=6.0, # chemotactic precision
    Γ=50.0, # chemotactic gain
    κ=50.0, # receptor gain (1/μM)
    τm=1.3, # chemotactic memory time (s)
    a=0.5, # bacterial radius (μm)
    rng=Xoshiro(1)
)
    space = ContinuousSpace(fill(L, SVector{3}); periodic)
    origin = fill(L/2, SVector{3})
    chemo = ExpField(origin, R, Cs, Cb, γ)
    source = Obstacles(space, [ObstacleSphere(origin, R)]; contact=HaltUntilTumble())
    properties = Dict(
        :chemoattractant => chemo,
        :obstacles => (source,),
    )
    model = StandardABM(microbe_type(mot), space, dt;
        properties,
        container=Vector,
        rng,
    )
    for i in 1:n
        pos = random_free_position(model; radius=a)
        add_agent!(pos, model;
            motility=make_motility(mot, U, λ),
            rotational_diffusivity=Drot,
            radius=a,
            behaviors=(Brumley(; chemotactic_precision=Π, gain=Γ, gain_receptor=κ, memory=τm),),
        )
    end
    return model
end

function setup_abm_community(;
    n=10_000, # number of bacteria
    L=1e3, # domain size (μm)
    periodic=true, # periodic boundary conditions (true / false)
    Rmin=0.65, # minimum phytoplankton radius in the community (μm)
    Rmax=10.0, # maximum phytoplankton radius in the community (μm)
    α=0.85, # phytoplankton allometric scaling
    Aphy=1e6, # total phytoplankton abundance (cells/mL)
    PER=0.5, # phytoplankton percent extracellular release
    Cb=0.0, # background chemoattractant (μM)
    γ=150.0, # chemoattractant decay lengthscale (μm)
    σ=40.0, # Ewald splitting length of the concentration field (μm)
    mot="RT", # motility pattern (RT / RR / RRF)
    dt=0.02, # integration timestep (s)
    U=45.0, # swimming speed (μm/s)
    λ=2.2, # reorientation rate (1/s)
    Drot=0.035, # rotational diffusivity (rad^2/s)
    Π=6.0, # chemotactic precision
    Γ=50.0, # chemotactic gain
    κ=50.0, # receptor gain (1/μM)
    τm=1.3, # chemotactic memory time (s)
    a=0.5, # bacterial radius (μm)
    teq=0.0, # equilibration time (s) before exposure starts accumulating
    rng=Xoshiro(1), # seed
)
    periodic || throw(ArgumentError("the Ewald community field requires a periodic domain"))
    space = ContinuousSpace(fill(L, SVector{3}); periodic)
    community = generate_community(space, α, Aphy, Rmin, Rmax, PER; rng)
    spatialcfg, phytoplankton_leakage = community
    phytoplankton_positions = SVector{3}.(getproperty.(spatialcfg, :pos))
    phytoplankton_radii = getproperty.(spatialcfg, :radius)
    system = ParticleSystem(
        xpositions=zeros(SVector{3,Float64}, n),
        ypositions=phytoplankton_positions,
        unitcell=spacesize(space),
        cutoff=3.5σ, # erfc(3.5) ≈ 7e-7
        output=OutCommField(zeros(n), zeros(SVector{3,Float64}, n)),
        output_name=:measurements,
        parallel=false
    )
    phytoplankton = Obstacles(space, ObstacleSphere.(phytoplankton_positions, phytoplankton_radii);
        contact=HaltUntilTumble()
    )
    chemo = CommField(Cb, γ, σ)
    strengths = point_strength.(phytoplankton_radii, phytoplankton_leakage, γ)
    properties = Dict(
        :neighborlist => system,
        :ewald_far => EwaldFar(phytoplankton_positions, strengths, L, γ, σ),
        :obstacles => (phytoplankton,),
        :phytoplankton_radii => phytoplankton_radii,
        :phytoplankton_leakage => phytoplankton_leakage,
        :chemoattractant => chemo,
        :exposure => zeros(n),
        :exposure_from => round(Int, teq/dt),
    )
    model = StandardABM(microbe_type(mot), space, dt;
        properties,
        container = Vector,
        agent_step! = dummystep,
        model_step! = community_step!,
        rng,
    )
    for i in 1:n
        pos = random_free_position(model; radius=a)
        add_agent!(pos, model;
            motility=make_motility(mot, U, λ),
            rotational_diffusivity=Drot,
            radius=a,
            behaviors=(Brumley(; chemotactic_precision=Π, gain=Γ, gain_receptor=κ, memory=τm),),
        )
    end
    community_field!(model)
    return model
end

setup_abm_community_Δche(; kwargs...) =
    setup_abm_community(κ=0.0, Γ=0.0; kwargs...)

using IterTools: partition
using StatsBase: geomean
using BubbleBath
function generate_community(space, α, Aphy, Rmin, Rmax, PER; rng=Xoshiro(1))
    Nclasses = 7
    sizeclasses = partition(logrange(Rmin, Rmax; length=Nclasses+1), 2, 1)
    # geometric mean of radii rounded to first decimal for convenience
    rs = round.(sigdigits=2, geomean.(sizeclasses))
    Nrel = map(r -> 1/r^(3α), rs) # relative abundance in each size class
    vol = prod(spacesize(space)) * 1e-12 # domain volume in mL
    Nperml = Nrel ./ sum(Nrel) .* Aphy # abundance per mL
    N = round.(Int, Nperml .* vol) # abundance in domain
    radii = reduce(vcat, [fill(r, n) for (r, n) in zip(rs, N)])
    spatialcfg = bubblebath(radii, Tuple(spacesize(space)); rng)
    intensities = leaked_concentration.(getproperty.(spatialcfg, :radius), PER)
    spatialcfg, intensities
end

using Unitful
function leakage_rate(R, PER;
    nC=5,
    μ=1.0*1u"d^-1",
    b=1.67e-4u"mol/cm^2.28",
    D=500u"μm^2/s"
)
    Ru = R * 1u"μm"
    L = b * Ru^2.28 * μ * PER / nC
    ustrip(L |> u"pmol/hr")
end
function leaked_concentration(R, PER;
    nC=5,
    μ=1.0*1u"d^-1",
    b=1.67e-4u"mol/cm^2.28",
    D=500u"μm^2/s"
)
    Ru = R * 1u"μm"
    L = b * Ru^2.28 * μ * PER / nC
    Cs = L / (4π*D*Ru) |> u"μM"
    ustrip(Cs)
end

function community_step!(model)
    # 1. move every swimmer (surface contacts resolved by the obstacles),
    #    then rotational diffusion
    for microbe in allagents(model)
        move_step!(microbe, model)
    end

    # 2. recompute concentration + gradient at the new positions
    community_field!(model)
    accumulate_exposure!(model)

    # 3. chemotactic update and reorientation
    for microbe in allagents(model)
        affect_step!(microbe, model)
        reorient_step!(microbe, model)
    end
end

# Ewald sum at every microbe position: near part over the neighbor list,
# then the far part from the precomputed mesh
function community_field!(model)
    nlist = abmproperties(model)[:neighborlist]
    for i in allids(model)
        nlist.xpositions[i] = position(model[i])
    end
    map_pairwise!(
        (x, y, i, j, d2, out) -> comm_out!(x, y, i, j, d2, out, model),
        nlist
    )
    far = abmproperties(model)[:ewald_far]
    out = nlist.measurements
    for i in allids(model)
        x = position(model[i])
        out.c[i] += far(x)
        out.dc[i] += far_gradient(far, x)
    end
    return model
end

# time-integrated excess concentration (μM·s) experienced by each microbe,
# summed over the states reached after the first `exposure_from` steps
function accumulate_exposure!(model)
    abmtime(model) < abmproperties(model)[:exposure_from] && return model
    exposure = abmproperties(model)[:exposure]
    c = abmproperties(model)[:neighborlist].measurements.c
    dt = abmtimestep(model)
    for i in allids(model)
        exposure[i] += dt * c[i]
    end
    return model
end

function make_motility(mot, U, λ)
    if mot == "RT"
        RunTumble([U], 1/λ, Isotropic(3))
    elseif mot == "RR"
        # keyword form: the positional RunReverse method of MicrobeAgents
        # does not forward keyword arguments
        RunReverse(;
            run_speed_forward=[U], run_duration_forward=1/λ,
            run_speed_backward=[U], run_duration_backward=1/λ,
            angle=Normal(-π, π/15),
        )
    elseif mot == "RRF"
        RunReverseFlick([U], 1/λ, [U], 1/λ;
            reverse_angle=Normal(-π, π/15),
            flick_angle=MixtureModel([Normal(π/2,π/6), Normal(-π/2,π/6)]),
        )
    else
        throw(ArgumentError("unknown motility pattern $mot"))
    end
end

# concrete agent type, so that the `Vector` container is concretely typed
microbe_type(mot) = Microbe{3, (mot == "RT" ? 2 : 4), Tuple{Brumley}}

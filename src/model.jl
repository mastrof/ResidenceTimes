export setup_abm, setup_abm_community
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
)
    space = ContinuousSpace(fill(L, SVector{3}); periodic)
    origin = fill(L/2, SVector{3})
    chemo = ExpField(origin, R, Cs, Cb, γ)
    properties = Dict(
        :chemoattractant => chemo
    )
    N = (mot == "RT" ? 2 : 4)
    model = StandardABM(Brumley{3,N}, space, dt;
        properties,
        container=Vector,
    )
    for i in 1:n
        motility = if mot == "RT"
            RunTumble([U], 1/λ, Isotropic(3))
        elseif mot == "RR"
            RunReverse([U], 1/λ, [U], 1/λ)
        elseif mot == "RRF"
            RunReverseFlick([U], 1/λ, [U], 1/λ)
        end
        add_agent!(model;
            motility,
            rotational_diffusivity=Drot,
            chemotactic_precision=Π,
            gain=Γ,
            gain_receptor=κ,
            memory=τm,
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
    mot="RT", # motility pattern (RT / RR / RRF)
    dt=0.02, # integration timestep (s)
    U=45.0, # swimming speed (μm/s)
    λ=2.2, # reorientation rate (1/s)
    Drot=0.035, # rotational diffusivity (rad^2/s)
    Π=6.0, # chemotactic precision
    Γ=50.0, # chemotactic gain
    κ=50.0, # receptor gain (1/μM)
    τm=1.3, # chemotactic memory time (s)
    rng=Xoshiro(1), # seed
)
    space = ContinuousSpace(fill(L, SVector{3}); periodic)
    community = generate_community(space, α, Aphy, Rmin, Rmax, PER; rng)
    spatialcfg, phytoplankton_leakage = community
    phytoplankton_positions = SVector{3}.(getproperty.(spatialcfg, :pos))
    phytoplankton_radii = getproperty.(spatialcfg, :radius)
    system = ParticleSystem(
        xpositions=zeros(SVector{3,Float64}, n),
        ypositions=phytoplankton_positions,
        unitcell=spacesize(space),
        cutoff=3γ,
        output=OutCommField(zeros(n), zeros(SVector{3,Float64}, n)),
        output_name=:measurements,
        parallel=false
    )
    chemo = CommField(Cb, γ)
    properties = Dict(
        :neighborlist => system,
        :phytoplankton_radii => phytoplankton_radii,
        :phytoplankton_leakage => phytoplankton_leakage,
        :chemoattractant => chemo,
    )
    N = (mot == "RT" ? 2 : 4)
    model = StandardABM(Brumley{3,N}, space, dt;
        properties,
        container = Vector,
        agent_step! = dummystep,
        model_step! = community_step!,
    )
    for i in 1:n
        motility = if mot == "RT"
            RunTumble([U], 1/λ, Isotropic(3))
        elseif mot == "RR"
            RunReverse([U], 1/λ, [U], 1/λ)
        elseif mot == "RRF"
            RunReverseFlick([U], 1/λ, [U], 1/λ)
        end
        add_agent!(model;
            motility,
            rotational_diffusivity=Drot,
            chemotactic_precision=Π,
            gain=Γ,
            gain_receptor=κ,
            memory=τm,
        )
    end
    abmproperties(model).neighborlist.xpositions .=
        position.(allagents(model))
    map_pairwise!(
        (x,y,i,j,d2,out) -> comm_out!(x,y,i,j,d2,out,model),
        abmproperties(model)[:neighborlist]
    )
    return model
end

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
function leaked_concentration(R, PER;
    nC=5,
    μ=1.0*1u"d^-1",
    b=1.67e-4u"mol/cm^2.28",
    D=500u"μm^2/s"
)
    Ru = R * 1u"μm"
    L = b * Ru^2.28 * μ * PER * nC
    Cs = L / (4π*D*Ru) |> u"μM"
    ustrip(Cs)
end

function community_step!(model)
    dt = abmtimestep(model)
    # move all cells to new positions
    for microbe in allagents(model)
        move_agent!(microbe, model, speed(microbe)*dt)
        MicrobeAgents.rotational_diffusion!(microbe, model)
    end
    # update neighbor list
    nlist = abmproperties(model)[:neighborlist]
    for i in allids(model)
        nlist.xpositions[i] = position(model[i])
    end
    # compute new values of concentration and gradient
    map_pairwise(
        (x,y,i,j,d2,out) -> comm_out!(x,y,i,j,d2,out,model),
        nlist
    )
    # continue with individual microbe steps as usual
    for microbe in allagents(model)
        model.affect!(microbe, model)
        if MicrobeAgents.can_turn(microbe)
            MicrobeAgents.turn!(microbe, model)
        end
        p = switching_probability(microbe, model)
        if rand(abmrng(model)) < p
            update_motilestate!(microbe, model)
            new_motilestate = motilestate(microbe)
            if variantof(new_motilestate) === TurnState && iszero(duration(new_motilestate))
                MicrobeAgents.turn!(microbe, model)
                update_motilestate!(microbe, model)
            end
            MicrobeAgents.update_speed!(microbe, model)
        end
    end
end

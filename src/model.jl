export setup_abm

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

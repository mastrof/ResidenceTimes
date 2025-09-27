module ResidenceTimes

using DrWatson
using LinearAlgebra
using Agents
using MicrobeAgents
using StaticArrays
using CSV
using DataFrames

include("utils.jl")
include("wsave.jl")
include("reader.jl")
include("concentration_field.jl")
include("model.jl")

end # module

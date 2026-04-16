module ResidenceTimes

using DrWatson
using LinearAlgebra
using Agents
using MicrobeAgents
using StaticArrays
using CellListMap
using CSV
using DataFrames
using Random

include("utils.jl")
include("wsave.jl")
include("reader.jl")
include("concentration_field.jl")
include("model.jl")

end # module

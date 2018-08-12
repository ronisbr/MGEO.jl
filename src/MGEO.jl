module MGEO

import Base: convert

using DelimitedFiles
using Distributed
using LinearAlgebra
using StaticArrays
using Random

export Design_Variable, MGEO_Structure, Pareto_Point

################################################################################
#                                    Files
################################################################################

include("types.jl")
include("configuration.jl")
include("misc.jl")
include("mgeo_canonical.jl")
include("mgeo_var.jl")

end # module

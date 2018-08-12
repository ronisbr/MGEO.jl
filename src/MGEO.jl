module MGEO

import Base: convert

using DelimitedFiles
using Distributed
using LinearAlgebra
using StaticArrays
using Random

################################################################################
#                                    Files
################################################################################

include("types.jl")
include("configuration.jl")
include("misc.jl")
include("mgeo_canonical.jl")
include("mgeo_var.jl")
include("mgeo_real.jl")

end # module

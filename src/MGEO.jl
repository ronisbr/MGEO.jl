module MGEO

import Base: convert

using DelimitedFiles
using Distributed
using LinearAlgebra
using StaticArrays
using Random

export Design_Variable, MGEO_Structure, Pareto_Point

################################################################################
#                                  Structures
################################################################################

#                              Public Structures
# ==============================================================================

"""
    struct Design_Variable{Nb}

Structure that defines the limits of the design variables.

# Fields

* `bits`: Number of bits.
* `min`: Minimum allowed value for the design variable.
* `max`: Maximum allowed value for the design variable.
* `full_scale`: Full scale of the design variable.
* `index`: Index in the string.
* `name`: Name of the variable.

"""
struct Design_Variable
    bits::Int64
    min::Float64
    max::Float64
    full_scale::Int64
    index::Int64
    name::AbstractString
end

"""
    struct Pareto_Point{Nv, Nf}

Structure that defines a point in the Pareto frontier.

`Nv` is the number of design variables and `Nf` is the number of objective
functions.

# Fields

* `vars`: Design variables.
* `f`: Objective functions.

"""
struct Pareto_Point{Nv, Nf}
    vars::SVector{Nv, Float64}
    f::SVector{Nf, Float64}
end

"""
    struct MGEO_Structure{Nv, Nf}

Structure to store the configuration of MGEO.

`Nv` is the number of design variables and `Nf` is the number of objective
functions.

# Fields

* `tau`: Parameter to set the search determinism.
* `ngen_max`: Maximum number of generations.
* `run_max`: Maximum number of independent runs (reinitializations).
* `design_vars`: Structure that store the configuration for each design variable.
* `mgeo_eps`: Epsilon to compare objective functions.

"""
struct MGEO_Structure{Nv, Nf}
    τ::Float64
    ngen_max::Int64
    run_max::Int64
    design_vars::SVector{Nv, Design_Variable}
    mgeo_eps::Float64
end

#                              Private Structures
# ==============================================================================

"""
    struct sRank

Structure used to sort the candidate points.

# Fields

* `valid`: Is the value valid?
* `index`: Index of the flipped bit.
* `f`: Objective function value after flipping the bit.

"""
struct sRank
    valid::Bool
    index::Int64
    f::Float64
end

################################################################################
#                                    Files
################################################################################

include("misc.jl")
include("run.jl")

end # module

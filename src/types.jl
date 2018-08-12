# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
#
# Description
#
#   Definition of the types and structures.
#
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #

export Design_Variable_MGEO_Canonical
export Design_Variable_MGEO_Var
export Design_Variable_MGEO_Real
export MGEO_Structure, Pareto_Point
export MGEO_Canonical, MGEO_Var, MGEO_Real

MGEO_Canonical() = Val{:MGEO_Canonical}
MGEO_Var()       = Val{:MGEO_Var}
MGEO_Real()      = Val{:MGEO_Real}

################################################################################
#                                  Structures
################################################################################

#                              Public Structures
# ==============================================================================

"""
    struct Design_Variable_MGEO_Canonical

Structure that defines the design variables for MGEO Canonical.

# Fields

* `bits`: Number of bits.
* `min`: Minimum allowed value for the design variable.
* `max`: Maximum allowed value for the design variable.
* `full_scale`: Full scale of the design variable.
* `index`: Index in the string.
* `name`: Name of the variable.

"""
struct Design_Variable_MGEO_Canonical
	bits::Int64
	min::Float64
	max::Float64
	full_scale::Int64
	index::Int64
	name::AbstractString
end

"""
    struct Design_Variable_MGEO_Var

Structure that defines the design variables for MGEO Var.

# Fields

* `bits`: Number of bits.
* `min`: Minimum allowed value for the design variable.
* `max`: Maximum allowed value for the design variable.
* `full_scale`: Full scale of the design variable.
* `index`: Index in the string.
* `name`: Name of the variable.

"""
struct Design_Variable_MGEO_Var
	bits::Int64
	min::Float64
	max::Float64
	full_scale::Int64
	index::Int64
	name::AbstractString
end

"""
    struct Design_Variable_MGEO_Real

Structure that defines the design variables for MGEO Real.

# Fields

* `min`: Minimum allowed value for the design variable.
* `max`: Maximum allowed value for the design variable.
* `σ`: Standard deviation used when modifying the variable.
* `name`: Name of the variable.

"""
struct Design_Variable_MGEO_Real
	min::Float64
	max::Float64
    σ::Float64
	name::AbstractString
end

"""
    struct MGEO_Structure{Nv, Nf, T}

Structure to store the configuration of MGEO.

`Nv` is the number of design variables, `Nf` is the number of objective
functions, and `T` is the MGEO algorithm that will be used.

# Fields

* `tau`: Parameter to set the search determinism.
* `ngen_max`: Maximum number of generations.
* `run_max`: Maximum number of independent runs (reinitializations).
* `design_vars`: Structure that store the configuration for each design variable.
* `mgeo_eps`: Epsilon to compare objective functions.

"""
struct MGEO_Structure{Nv, Nf, T}
    τ::Float64
    ngen_max::Int64
    run_max::Int64
    design_vars::SVector{Nv, T}
    mgeo_eps::Float64
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

#                              Private Structures
# ==============================================================================

"""
    struct sRank

Structure used to sort the candidate points in MGEO Canonical and MGEO Var.

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

"""
    struct sRank_Real

Structure used to sort the candidate points in MGEO Real.

# Fields

* `valid`: Is the value valid?
* `index`: Index of the modified variable.
* `var_value`: Value of the modified variable.
* `f`: Objective function value after flipping the bit.

"""
struct sRank_Real
    valid::Bool
    index::Int64
    var_value::Float64
    f::Float64
end


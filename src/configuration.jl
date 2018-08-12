# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
#
# Description
#
#   Functions to configure the MGEO.
#
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #

################################################################################
#                              Public Functions
################################################################################

"""
    function conf_design_vars(T, Nv::Int64, bits::Int64, min::Number, max::Number)

This function configures `Nv` design variables in which all of them will have
the same number of bits `bits`, and maximum and minimum values `min` and `max`.

The parameter `T` selects which MGEO algorithm will be used and it can be
`MGEO_Canonical()` (**Default**) or `MGEO_Var()`.

    function conf_design_vars(Val{:MGEO_Real}, Nv::Int64, min::Number, max::Number, σ::Number)

This function configures `Nv` design variables in which all of them will have
the same maximum and minimum values `min` and `max`, and the same standard
deviation `σ` used when modifying the variables. This is used for MGEO Real
algorithm.

# Args

* `T`: (OPTIONAL) Type of algorithm (**Default**: `MGEO_Canonical()`).
* `Nv`: The number of design variables.
* `bits`: The number of bits for each variable.
* `min`: The minimum value for all the design variables.
* `max`: The maximum value for all the design variables.

* `σ`: Standard deviation used when modifying a variable in MGEO Real.

# Returns

An array with the configured design variables.

"""
conf_design_vars(Nv::Int64, bits::Int64, min::Number, max::Number) =
    conf_design_vars(Val{:MGEO_Canonical}, Nv, bits, min, max)

function conf_design_vars(T::Type{Val{:MGEO_Canonical}},
                          Nv::Int64,
                          bits::Int64,
                          min::Number,
                          max::Number)

    # Check if the number of design variables is greater than 1.
    ( Nv <= 0 ) && @error("The number of design variables must be higher than 0.")

    # Check if the number of bits are greater than 0.
    ( bits <= 0 ) && @error("The number of bits must be higher than 0.")

    # Check if the minimum value is smaller than maximum value.
    ( min >= max ) && @error("The minimum value must be less than the maximum value for a variable.")

    # Create the array of design variables.
    design_vars = Vector{Design_Variable_MGEO_Canonical}(undef, Nv)

    # Full scale for each variable.
    fullscale = (1 << bits) - 1

    @inbounds for i=1:Nv
        design_vars[i] = Design_Variable_MGEO_Canonical(bits,
                                                        min,
                                                        max,
                                                        fullscale,
                                                        (i-1)*bits+1,
                                                        "Var. $i")
    end

    SVector{Nv, Design_Variable_MGEO_Canonical}(design_vars)
end

function conf_design_vars(T::Type{Val{:MGEO_Var}},
                          Nv::Int64,
                          bits::Int64,
                          min::Number,
                          max::Number)

    # Check if the number of design variables is greater than 1.
    ( Nv <= 0 ) && @error("The number of design variables must be higher than 0.")

    # Check if the number of bits are greater than 0.
    ( bits <= 0 ) && @error("The number of bits must be higher than 0.")

    # Check if the minimum value is smaller than maximum value.
    ( min >= max ) && @error("The minimum value must be less than the maximum value for a variable.")

    # Create the array of design variables.
    design_vars = Vector{Design_Variable_MGEO_Var}(undef, Nv)

    # Full scale for each variable.
    fullscale = (1 << bits) - 1

    @inbounds for i=1:Nv
        design_vars[i] = Design_Variable_MGEO_Var(bits,
                                                  min,
                                                  max,
                                                  fullscale,
                                                  (i-1)*bits+1,
                                                  "Var. $i")
    end

    SVector{Nv, Design_Variable_MGEO_Var}(design_vars)
end

function conf_design_vars(T::Type{Val{:MGEO_Real}},
                          Nv::Int64,
                          min::Number,
                          max::Number,
                          σ::Number)

    # Check if the number of design variables is greater than 1.
    ( Nv <= 0 ) && @error("The number of design variables must be higher than 0.")

    # Check if the minimum value is smaller than maximum value.
    ( min >= max ) && @error("The minimum value must be less than the maximum value for a variable.")

    # Check if standard deviation σ is higher than 0.
    ( σ <= 0 ) && @error("The standard deviation must be higher than 0.")

    # Create the array of design variables.
    design_vars = Vector{Design_Variable_MGEO_Real}(undef, Nv)

    @inbounds for i=1:Nv
        design_vars[i] = Design_Variable_MGEO_Real(min, max, σ, "Var. $i")
    end

    SVector{Nv, Design_Variable_MGEO_Real}(design_vars)
end

"""
    function conf_design_vars(T, bits::Vector{Int64}, min::Vector{T1}, max::Vector{T2}, var_names::Vector{String}) where {T1<:Number, T2<:Number}

This function configures the design variables specifying for each one the number
of bits `bits`, the minimum `min` and maximum `max` values, and the variable
names `var_names`.

The parameter `T` selects which MGEO algorithm will be used and it can be
`MGEO_Canonical()` (**Default**) or `MGEO_Var()`.

    function conf_design_vars(T::Type{Val{:MGEO_Real}}, min::Vector{T1}, max::Vector{T2}, σ::Vector{T3}, var_names::Vector{String}) where

This function configures the design variables specifying for each one the
minimum `min` and maximum `max` values, the standard deviations `σ` used when
modifying the variable, and the variable names `var_names`.

# Args

* `T`: (OPTIONAL) Type of algorithm (**Default**: `MGEO_Canonical()`).
* `bits`: List containing the number of bits for each design variable.
* `min`: List containing the minimum of each design variable.
* `max`: List containing the maximum of each design variable.
* `var_names`: List containing the name of each design variable.

* `σ`: List containing the standard deviations used when modifying a variable in
       MGEO Real.

# Returns

An array with the configured design variables.

"""
conf_design_vars(bits::Vector{Int64},
                 min::Vector{T1},
                 max::Vector{T2},
                 var_names::Vector{String}) where
    {T1<:Number, T2<:Number} =
    conf_design_vars(Val{:MGEO_Canonical}, bits, min, max, var_names)

function conf_design_vars(T::Type{Val{:MGEO_Canonical}},
                          bits::Vector{Int64},
                          min::Vector{T1},
                          max::Vector{T2},
                          var_names::Vector{String}) where
 	{T1<:Number, T2<:Number}

    # Check if the size of arrays is correct.
    Nv = length(bits)

    ( ( length(min)       != Nv ) ||
      ( length(max)       != Nv ) ||
      ( length(var_names) != Nv ) ) &&
    @error("The size of the vectors does not match.")

    # Number of bits.
    num_bits = 0

    # Create the array of design variables.
    design_vars = Vector{Design_Variable{T}}(undef, Nv)

    for i = 1:Nv
        # Check if the number of bits are larger than 0.
        ( bits[i] <= 0 ) && @error("The number of bits must be higher than 0.")

        # Check if the min value is smaller than max.
        ( min[i] >= max[i] ) && @error("The minimum value must be less than the maximum value for a variable.")

        design_vars[i] = Design_Variable_MGEO_Canonical(bits[i],
                                                        min[i],
                                                        max[i],
                                                        (UInt64(1) << bits[i]) - 1,
                                                        num_bits+1,
                                                        var_names[i])
        num_bits += bits[i]
    end

    SVector{Nv, Design_Variable_MGEO_Canonical}(design_vars)
end

function conf_design_vars(T::Type{Val{:MGEO_Var}},
                          bits::Vector{Int64},
                          min::Vector{T1},
                          max::Vector{T2},
                          var_names::Vector{String}) where
 	{T1<:Number, T2<:Number}

    # Check if the size of arrays is correct.
    Nv = length(bits)

    ( ( length(min)       != Nv ) ||
      ( length(max)       != Nv ) ||
      ( length(var_names) != Nv ) ) &&
    @error("The size of the vectors does not match.")

    # Number of bits.
    num_bits = 0

    # Create the array of design variables.
    design_vars = Vector{Design_Variable{T}}(undef, Nv)

    for i = 1:Nv
        # Check if the number of bits are larger than 0.
        ( bits[i] <= 0 ) && @error("The number of bits must be higher than 0.")

        # Check if the min value is smaller than max.
        ( min[i] >= max[i] ) && @error("The minimum value must be less than the maximum value for a variable.")

        design_vars[i] = Design_Variable_MGEO_Var(bits[i],
                                                  min[i],
                                                  max[i],
                                                  (UInt64(1) << bits[i]) - 1,
                                                  num_bits+1,
                                                  var_names[i])
        num_bits += bits[i]
    end

    SVector{Nv, Design_Variable_MGEO_Var}(design_vars)
end

function conf_design_vars(T::Type{Val{:MGEO_Real}},
                          min::Vector{T1},
                          max::Vector{T2},
                          σ::Vector{T3},
                          var_names::Vector{String}) where
 	{T1<:Number, T2<:Number, T3<:Number}

    # Check if the size of arrays is correct.
    Nv = length(min)

    ( ( length(min)       != Nv ) ||
      ( length(max)       != Nv ) ||
      ( length(σ)         != Nv ) ||
      ( length(var_names) != Nv ) ) &&
    @error("The size of the vectors does not match.")

    # Create the array of design variables.
    design_vars = Vector{Design_Variable{T}}(undef, Nv)

    for i = 1:Nv
        # Check if the min value is smaller than max.
        ( min[i] >= max[i] ) && @error("The minimum value must be less than the maximum value for a variable.")

        # Check if standard deviation σ is higher than 0.
        ( σ[i] <= 0 ) && @error("The standard deviation must be higher than 0.")

        design_vars[i] = Design_Variable_MGEO_Real(min[i],
                                                   max[i],
                                                   σ[i],
                                                   var_names[i])
    end

    SVector{Nv, Design_Variable_MGEO_Real}(design_vars)
end

"""
    function conf_design_vars(bits::Vector{Int64}, min::Vector{T1}, max::Vector{T2}) where {T1<:Number, T2<:Number}

This function configures the design variables specifying for each one the number
of bits `bits`, and the minimum `min` and maximum `max` values. The variable
names will be selected as `Var. 1`, `Var. 2`, etc.

The parameter `T` selects which MGEO algorithm will be used and it can be
`MGEO_Canonical()` (**Default**) or `MGEO_Var()`.

    function conf_design_vars(T::Type{Val{:MGEO_Real}}, min::Vector{T1}, max::Vector{T2}, σ::Vector{T3}) where {T1<:Number, T2<:Number, T3<:Number}

This function configures the design variables specifying for each one the
minimum `min` and maximum `max` values, and the standard deviations `σ` used
when modifying the variable. The variable names will be selected as `Var. 1`,
`Var. 2`, etc.

# Args

* `T`: (OPTIONAL) Type of algorithm (**Default**: `MGEO_Canonical()`).
* `bits`: List containing the number of bits for each design variable.
* `min`: List containing the minimum of each design variable.
* `max`: List containing the maximum of each design variable.

* `σ`: List containing the standard deviations used when modifying a variable in
       MGEO Real.

# Returns

An array of `Design_Variable` with the specified values.

"""
conf_design_vars(bits::Vector{Int64},
                 min::Vector{T1},
                 max::Vector{T2}) where {T1<:Number, T2<:Number} =
    conf_design_vars(Val{:MGEO_Canonical}, bits, min, max)

function conf_design_vars(T::Union{Type{Val{:MGEO_Canonical}},
								   Type{Val{:MGEO_Var}}},
                          bits::Vector{Int64},
                          min::Vector{T1},
                          max::Vector{T2}) where {T1<:Number, T2<:Number}

    # Create an array with variable names.
    var_names = Array(String, length(bits))

    for i=1:length(varNames)
        varNames[i] = "Var. $i"
    end

    conf_design_vars(T, bits, min, max, var_names)
end

function conf_design_vars(T::Type{Val{:MGEO_Real}},
                          min::Vector{T1},
                          max::Vector{T2},
                          σ::Vector{T3}) where
    {T1<:Number, T2<:Number, T3<:Number}

    # Create an array with variable names.
    var_names = Array(String, length(bits))

    for i=1:length(varNames)
        varNames[i] = "Var. $i"
    end

    conf_design_vars(T, min, max, σ, var_names)
end

"""
    function conf_mgeo(Nf::Int64, τ::Float64, ngen_max::Int64, run_max::Int64, design_vars::SVector{Nv, Design_Variable{T}}, mgeo_eps::Float64 = 1e-10) where {Nv,T}

Configure the MGEO.

The MGEO algorithm that will be used will be the same configured in
`design_vars`.

# Args

* `Nf`: Number of objective functions.
* `τ`: Parameter τ of MGEO.
* `ngen_max`: Maximum number of generations.
* `run_max`: Maximum number of independent runs.
* `design_vars`: Configuration of the design variables (see `conf_design_vars`).
* `mgeo_eps`: (OPTIONAL) Epsilon to compare objective functions (**Default**:
              `1e-10`)

# Returns

An instance of the structure `MGEO_Structure` with the MGEO configuration.

"""
function conf_mgeo(Nf::Int64,
                   τ::Float64,
                   ngen_max::Int64,
                   run_max::Int64,
                   design_vars::SVector{Nv, T},
                   mgeo_eps::Float64 = 1e-10) where {Nv,T}

    # Check if the number of objective functions is higher than 1.
    ( Nf <= 0 ) && @error("The number of objective functions must be higher than 0.")

    MGEO_Structure{Nv, Nf, T}(τ, ngen_max, run_max, design_vars, mgeo_eps)
end

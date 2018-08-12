# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
#
# Description
#
#   Auxiliary functions for MGEO.
#
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #

export conf_design_vars, conf_mgeo, convert, pf_to_csv, sort_pareto!

################################################################################
#                                  Overloads
################################################################################

"""
    function convert(::Type{Matrix}, pareto::Vector{Pareto_Point{Nv,Nf}}) where {Nv,Nf}

Convert the structure that stores the Pareto frontier `pareto` to an array.

# Args

* `paretoFrontier`: Pareto frontier.

# Returns

A matrix of `Float64` in which the i-th line of the array is

    [var[1] var[2] ... var[n] f[1] f[2] ... f[nf]],

where `n` is the number of design variables and `nf` is the number of
objective functions.

"""
function convert(::Type{Matrix},
                 pareto::Vector{Pareto_Point{Nv,Nf}}) where {Nv,Nf}

    # Get the number of points.
    np = length(pareto)

    # Allocate the array.
    pf_array = zeros(np, Nv+Nf)

    # Fill the array.
    for i=1:np
        pf_array[i,:] = [pareto[i].vars' pareto[i].f']
    end

    pf_array
end

################################################################################
#                              Public Functions
################################################################################

"""
    function conf_design_vars(Nv::Int64, bits::Int64, min::Number, max::Number)

Configure `Nv` design variables in which all of them will have the same number
of bits `bits`, and maximum and minimum values `min` and `max`.

# Args

* `Nv`: The number of design variables.
* `bits`: The number of bits for each variable.
* `min`: The minimum value for all the design variables.
* `max`: The maximum value for all the design variables.

# Returns

An array of `Design_Variable` with the specified values.

"""
function conf_design_vars(Nv::Int64, bits::Int64, min::Number, max::Number)
    # Check if the number of design variables is greater than 1.
    ( Nv <= 0 ) && @error("The number of design variables must be higher than 0.")

    # Check if the number of bits are greater than 0.
    ( bits <= 0 ) && @error("The number of bits must be higher than 0.")

    # Check if the minimum value is smaller than maximum value.
    ( min >= max ) && @error("The minimum value must be less than the maximum value for a variable.")

    # Create the array of design variables.
    design_vars = Vector{Design_Variable}(undef, Nv)

    # Full scale for each variable.
    fullscale = (1 << bits) - 1

    @inbounds for i=1:Nv
        design_vars[i] = Design_Variable(bits,
                                         min,
                                         max,
                                         fullscale,
                                         (i-1)*bits+1,
                                         "Var. $i")
    end

    SVector{Nv, Design_Variable}(design_vars)
end

"""
    function conf_design_vars(bits::Vector{Int64}, min::Vector{T1}, max::Vector{T2}, var_names::Vector{String}) where {T1<:Number, T2<:Number}

Configure the design variables.

# Args

* `bits`: List containing the number of bits for each design variable.
* `min`: List containing the minimum of each design variable.
* `max`: List containing the maximum of each design variable.
* `var_names`: List containing the name of each design variable.

# Returns

An array of `Design_Variable` with the specified values.

"""
function conf_design_vars(bits::Vector{Int64},
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
    design_vars = Vector{Design_Variable}(undef, Nv)

    for i = 1:Nv
        # Check if the number of bits are larger than 0.
        ( bits[i] <= 0 ) && @error("The number of bits must be higher than 0.")

        # Check if the min value is smaller than max.
        ( min[i] >= max[i] ) && @error("The minimum value must be less than the maximum value for a variable.")

        design_vars[i] = Design_Variable(bits[i],
                                         min[i],
                                         max[i],
                                         (UInt64(1) << bits[i]) - 1,
                                         num_bits+1,
                                         var_names[i])
        num_bits += bits[i]
    end

     design_vars
end

"""
    function conf_design_vars(bits::Vector{Int64}, min::Vector{T1}, max::Vector{T2}) where {T1<:Number, T2<:Number}

Configure the design variables.

# Args

* `bits`: List containing the number of bits for each design variable.
* `min`: List containing the minimum of each design variable.
* `max`: List containing the maximum of each design variable.

# Returns

An array of `Design_Variable` with the specified values.

"""
function conf_design_vars(bits::Vector{Int64},
                          min::Vector{T1},
                          max::Vector{T2}) where {T1<:Number, T2<:Number}

    # Create an array with variable names.
    var_names = Array(String, length(bits))

    for i=1:length(varNames)
        varNames[i] = "Var. $i"
    end

    conf_design_vars(bits, min, max, var_names)
end

"""
    function conf_mgeo(Nf::Int64, τ::Float64, ngen_max::Int64, run_max::Int64, design_vars::SVector{Nv, Design_Variable}, mgeo_eps::Float64 = 1e-10) where {Nv}

Configure the MGEO.

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
                   design_vars::SVector{Nv, Design_Variable},
                   mgeo_eps::Float64 = 1e-10) where {Nv}

    # Check if the number of objective functions is higher than 1.
    ( Nf <= 0 ) && @error("The number of objective functions must be higher than 0.")

    MGEO_Structure{Nv, Nf}(τ, ngen_max, run_max, design_vars, mgeo_eps)
end

"""
    function pf_to_csv(pareto::Vector{Pareto_Point{Nv,Nf}}, filename::String, design_vars::Vector{Design_Variable}) where {Nv,Nf}

Save the Pareto frontier `pareto` to a CSV file `filename`.

# Args

* `pareto`: Pareto frontier (array of `Pareto_Point`).
* `filename`: File name.
* `design_vars`: Information about design variables.

"""
function pf_to_csv(pareto::Vector{Pareto_Point{Nv,Nf}},
                   filename::String,
                   design_vars::SVector{Nv,Design_Variable}) where {Nv,Nf}

    # Convert Pareto frontier.
    pf_array = convert(Matrix, pareto)

    # Array of variable names.
    var_names = Matrix{String}(undef, 1, Nv)

    for i = 1:Nv
        var_names[1,i] = design_vars[i].name
    end

    # Array of function names.
    func_names = Matrix{String}(undef, 1, Nf)

    for i = 1:Nf
        func_names[1,i] = "Obj. Func. $i"
    end

    # Write the information to the file.
    writedlm(filename, [ [ var_names func_names ]; pf_array ], ",")

    nothing
end

"""
    function sort_pareto!(pareto::Vector{Pareto_Point{Nv,Nf}}, fobj::Int64) where {Nv,Nf}

Sort the points in Pareto frontier `pareto` using the objective function `fobj`.

# Args

* `pareto`: Pareto frontier.
* `fobj`: Objective function used to sort.

"""
function sort_pareto!(pareto::Vector{Pareto_Point{Nv,Nf}},
                      fobj::Int64) where {Nv,Nf}

    # Check the  input.
    nf = length(pareto[1].f)

    ( (fobj < 1) || (fobj > nf) ) && @error("fobj must be between 1 and the number of objective functions.")

    sort!(pareto, lt=(a,b)->begin
              # Check if data is valid.
              return (a.f[fobj] < b.f[fobj])
          end
         )

    nothing
end

################################################################################
#                              Private Functions
################################################################################

"""
    function bitstrtonum(design_vars::SVector{Nv, Design_Variable}, string::BitArray, factors::AbstractVector) where Nv

Convert the bit string `string` to real values using the information of the
design variables in `design_vars`.

# Args

* `design_vars`: Array with the description of design variables.
* `string`: Bit string.
* `factors`: Integer factors to convert the variable.

# Returns

An array of `Float64` with the real values of the design variables.

"""
function bitstrtonum(design_vars::SVector{Nv, Design_Variable},
                     string::BitArray,
                     factors::AbstractVector) where Nv

    # Array of variables.
    vars = Vector{Float64}(undef, Nv)

    # Loop for each design variable.
    @inbounds for i = 1:Nv
        # Convert the bits to integer.
        string_index_i = design_vars[i].index
        string_index_f = design_vars[i].index+design_vars[i].bits-1
        var_int = dot(factors[1:design_vars[i].bits],
                      string[string_index_i:string_index_f])

        # Convert the integer to real number.
        Δ       = (design_vars[i].max - design_vars[i].min)*var_int
        vars[i] = design_vars[i].min + Δ/(design_vars[i].full_scale)
    end

    SVector{Nv, Float64}(vars)
end

"""
    function check_dominance!(pareto::Vector{Pareto_Point{Nv,Nf}}, candidate::Pareto_Point{Nv,Nf}, mgeod::MGEO_Structure{Nv,Nf}) where {Nv, Nf}

Check the dominance of the point `candidate` in the Pareto frontier `pareto`.

# Args

* `pareto`: Pareto frontier.
* `candidate`: Candidate point to be added in the frontier.
* `mgeod`: MGEO structure (see `MGEO_Structure`).

# Returns

`true` if the point was added, `false` otherwise.

"""
function check_dominance!(pareto::Vector{Pareto_Point{Nv,Nf}},
                          candidate::Pareto_Point{Nv,Nf},
                          mgeod::MGEO_Structure{Nv,Nf}) where {Nv, Nf}
    # If Pareto frontier is empty, then add the point and return true.
    if length(pareto) == 0
        push!(pareto, candidate)
        return true
    end

    remove = Int64[]

    # Loop through the entire list of Pareto points.
    add_point = true

    for p=1:length(pareto)
        # Variable to check if the point in the Pareto frontier list is
        # dominated by the candidate.
        point_dominated     = true
        candidate_dominated = true

        # Loop through objective functions.
        @inbounds for i=1:Nf
            # If both objective functions are equal (given `mgeo_eps`), then
            # nothing can be stated about the dominance of both points.
            if abs(candidate.f[i] - pareto[p].f[i]) < mgeod.mgeo_eps
                continue
            end

            if candidate.f[i] < pareto[p].f[i]
                candidate_dominated = false
            elseif candidate.f[i] > pareto[p].f[i]
                point_dominated = false
            end
        end

        # If the point is dominated by the candidate and dominates the candidate
        # at the same time, then the candidate and the point is the same. Thus,
        # mark to do not add the candidate and exit the loop.

        if (point_dominated && candidate_dominated)
            add_point = false
            break
        end

        # If the point in the Pareto frontier is dominated, then mark it to be
        # excluded.
        if point_dominated
            push!(remove, p)
        # If the candidate is dominated by any point in the list, stop the
        # search and do not add it.
        elseif candidate_dominated
            add_point = false
            break
        end
    end

    # Remove the points dominated by the candidate point.
    @inbounds deleteat!(pareto, remove)

    # Check if the point must be added.
    add_point && push!(pareto, candidate)

    # Return if a point was added.
    add_point
end


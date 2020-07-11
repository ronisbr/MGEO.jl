# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
#
# Description
# ==============================================================================
#
#   Auxiliary functions for MGEO.
#
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #

export conf_design_vars, conf_mgeo, convert, pf_to_csv, sort_pareto!

################################################################################
#                                  Overloads
################################################################################

"""
    convert(::Type{Matrix}, pareto::Vector{Pareto_Point{Nv,Nf}}) where {Nv,Nf}

Convert the structure that stores the Pareto frontier `pareto` to an array.

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

"""
    pf_to_csv(pareto::Vector{Pareto_Point{Nv,Nf}}, filename::String, design_vars::SVector{Nv,T}) where {Nv,Nf,T}

Save the Pareto frontier `pareto` to a CSV file `filename` considering the
configuration of the design variables in `design_vars`.

"""
function pf_to_csv(pareto::Vector{Pareto_Point{Nv,Nf}},
                   filename::String,
                   design_vars::SVector{Nv,T}) where {Nv,Nf,T}

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
    sort_pareto!(pareto::Vector{Pareto_Point{Nv,Nf}}, fobj::Int64) where {Nv,Nf}

Sort the points in Pareto frontier `pareto` using the objective function index
`fobj`. Notice that the sort will be performed in place, modifying the vector
`pareto`.

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
    bitstrtonum(design_vars::SVector{Nv, T}, string::BitArray, factors::AbstractVector) where {Nv,T}

Convert the bit string `string` to real values using the information of the
design variables in `design_vars` and the integer factors `factors`.

# Returns

An array of `Float64` with the real values of the design variables.

"""
function bitstrtonum(design_vars::SVector{Nv, T},
                     string::BitArray,
                     factors::AbstractVector) where {Nv,T}

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
    check_dominance!(pareto::Vector{Pareto_Point{Nv,Nf}}, candidate::Pareto_Point{Nv,Nf}, mgeod::MGEO_Structure{Nv,Nf}) where {Nv, Nf}

Check the dominance of the point `candidate` in the Pareto frontier `pareto`
using the MGEO configuration in `mgeod`. Notice that the variable `pareto` can
be modified depending on the dominance of the `candidate`.

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

    @inbounds for p=1:length(pareto)
        # Variable to check if the point in the Pareto frontier list is
        # dominated by the candidate.
        point_dominated     = true
        candidate_dominated = true

        # Loop through objective functions.
        for i=1:Nf
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

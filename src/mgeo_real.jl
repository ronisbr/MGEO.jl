# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
#
# Description
#
#   Function to run the MGEO Real.
#
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #

export mgeo_run

"""
    function mgeo_run(mgeod::MGEO_Structure{Nv, Nf, Design_Variable_MGEO_Real}, f_obj::Function, show_debug::Bool = false) where {Nv, Nf}

Run the MGEO Real configured in `mgeod` using the objective functions `f_obj`.

# Args

* `mgeod`: Structure with the configurations of MGEO (see `MGEO_Structure`).
* `f_obj`: Objective functions.

# Keywords

* `show_debug`: Print debug information (**Default** = `true`).

# Returns

The Pareto frontier, which is an array of `Pareto_Point`.

##### Remarks

The objective function must have the following signature

        function f_obj(vars)
            1) Compute f[1], f[2], ..., f[nf] using vars[1], vars[2], ... var[n]
            2) return (valid_point, f)
        end

where `nf` is the number of objective functions, `n` is the number of design
variables, and `valid_point` is a boolean value that indicates if vars yield to
a valid point.

"""
function mgeo_run(mgeod::MGEO_Structure{Nv, Nf, Design_Variable_MGEO_Real},
                  f_obj::Function,
                  show_debug::Bool = false) where {Nv, Nf}
    # Pareto frontier.
    pareto = Pareto_Point{Nv,Nf}[]

    # Maximum number of generations per run.
    ngen_max_per_run = floor(mgeod.ngen_max/mgeod.run_max)

    # Vector to store the population.
    vars = Vector{Float64}(undef, Nv)

    # Initialization of Pareto Frontier
    # =================================

    # Initialize the string with one valid point in the Pareto Frontier.
    @inbounds for run = 1:mgeod.run_max
        # Sample a new population using a uniform distribution between the
        # minimum and maximum allowed range for each design variable.
        vars = Vector(map(x->x.min + rand()*(x.max-x.min), mgeod.design_vars))

        # Call the objective functions.
        (valid, f) = f_obj(vars)

        # Check if the number of objective functions are valid.
        (length(f) != Nf) && @error("The number of objective functions returned by f_obj is not equal to the one specified in mgeod.")

        if valid
            # Add the point to the Pareto Frontier.
            push!(pareto, Pareto_Point{Nv, Nf}(vars, SVector{Nf, Float64}(f)))
            break
        end
    end

    # Check if a valid point was found.
    if (length(pareto) == 0)
        return pareto
    end

    # Loop - Independent Runs
    # =======================

    @inbounds for run = 1:mgeod.run_max
        if show_debug
            println("--------------------------------------------------------------------------------")
            println("RUN = $run")
            println("")
            println("Generations ($ngen_max_per_run):")
            print("    0%")
        end

        # Sample a new string if it is not the first run.
        if run > 1
            vars = Vector(map(x->x.min + rand()*(x.max-x.min), mgeod.design_vars))
        end

        # Loop - MGEO Generations
        # =======================

        genDebug = 0

        @inbounds for ngen_per_run = 1:ngen_max_per_run
            if (show_debug)
                genDebug += 1

                if (genDebug >= 0.1*ngen_max_per_run)
                    genDebug = 0
                    print("...$(ngen_per_run/ngen_max_per_run*100)%")
                end
            end

            # Choose which objective function will be used to compute the
            # adaptability and to assemble the rank.
            chosen_func = rand(1:Nf)

            # List of all points created after the modification of population.
            candidate_points = Vector{Pareto_Point{Nv,Nf}}(undef,Nv)

            # Array to sort the candidate points.
            f_rank = Vector{sRank_Real}(undef, Nv)

            # Evaluate the objective functions using parallel computing.
            for i=1:Nv
                # Modify the i-th design variable.
                vars_i               = copy(vars)
                @inbounds vars_i[i] += mgeod.design_vars[i].σ*randn()

                # Evaluate the objective functions.
                (valid, f) = f_obj(vars_i)

                # Check if the solution is valid.
                if valid
                    # Create the candidate point.
                    candidate_point =
                        Pareto_Point{Nv, Nf}(vars, SVector{Nf, Float64}(f))

                    # Add the result to the rank.
                    @inbounds f_rank[i] = sRank_Real(true,
                                                     i,
                                                     vars_i[i],
                                                     f[chosen_func])
                    @inbounds candidate_points[i] = candidate_point
                else
                    @inbounds f_rank[i] = sRank_Real(false, i, vars_i[i], 0.0)
                end
            end

            # Add the points to the Pareto frontier.
            num_valid_points = 0

            @inbounds for i=1:length(candidate_points)
                if f_rank[i].valid
                    num_valid_points += 1
                    check_dominance!(pareto, candidate_points[i], mgeod)
                end
            end

            # If there are no valid points, an independent run must be called.
            (num_valid_points == 0) && break

            # Ranking
            # ========
            #
            # The points will be ranked according to the selected objective
            # function (`chosen_func`).
            #
            # Notice that the invalid points will be placed after the valid
            # ones.

            sort!(f_rank, lt=(a,b)->begin
                      # Check if data is valid.
                      if (a.valid && b.valid)
                          return a.f < b.f
                      elseif (a.valid)
                          return true
                      else
                          return false
                      end
                  end
                 )

            # Choose a variable to be changed for the new generation.
            var_accepted = false

            while !var_accepted
                # Sample a valid point.
                v_sample = rand(1:num_valid_points)

                # Accept the change with probability r^(-τ), where r is the
                # rank of the sampled variable.
                Pk = v_sample^(-mgeod.τ)

                if rand() <= Pk
                    # The bit is accepted, then exit the loop.
                    @inbounds vars[f_rank[v_sample].index] =
                        f_rank[v_sample].var_value
                    var_accepted = true
                end
            end
        end

        if show_debug
            println("")
            println("")
            println("Pareto frontier information:")
            println("    Number of points: $(length(pareto))")
            println("")
        end
    end

    pareto
end

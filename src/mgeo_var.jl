# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
#
# Description
#
#   Function to run the MGEO Var.
#
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #

export mgeo_run

"""
    function mgeo_run(mgeod::MGEO_Structure{Nv, Nf, Val{:MGEO_Var}}, f_obj::Function, show_debug::Bool = false) where {Nv, Nf}

Run the MGEO Var configured in `mgeod` using the objective functions `f_obj`.

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
function mgeo_run(mgeod::MGEO_Structure{Nv, Nf, Val{:MGEO_Var}},
                  f_obj::Function,
                  show_debug::Bool = false) where {Nv, Nf}
    # Pareto frontier.
    pareto = Pareto_Point{Nv,Nf}[]

    # Maximum number of generations per run.
    ngen_max_per_run = floor(mgeod.ngen_max/mgeod.run_max)

    # Get the number of bits the string must have by summing the bits for each
    # design variable. Also, get the maximum number of bits in one variable so
    # that we can compute the `factors` vector to convert the string of bits to
    # a real number.
    num_bits    = 0
    max_bit_num = 0

    @inbounds for k = 1:Nv
        num_bits += mgeod.design_vars[k].bits
        max_bit_num = max(max_bit_num, mgeod.design_vars[k].bits)
    end

    # Compute the vector with the interger factor to convert the string of bits
    # into a real number.
    factors = SVector{max_bit_num, Int64}(2 .^collect(0:1:(max_bit_num-1)))

    # Create the string.
    string = BitArray(undef,num_bits)

    # Initialization of Pareto Frontier
    # =================================

    # Initialize the string with one valid point in the Pareto Frontier.
    @inbounds for run = 1:mgeod.run_max
        # Sample a new string.
        rand!(string)

        # Convert string to real numbers.
        vars = bitstrtonum(mgeod.design_vars, string, factors)

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
        (run > 1) && rand!(string)

        # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
        #                   Loop - MGEO Generations
        # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #

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

            # List of all points created after flipping the bits.
            candidate_points = Vector{Pareto_Point{Nv,Nf}}(undef,num_bits)

            # Array to sort the candidate points.
            f_rank = Vector{sRank}(undef, num_bits)

            # Evaluate the objective functions using parallel computing.
            for i=1:num_bits
                # Flip the i-th bit.
                string_i              = copy(string)
                @inbounds string_i[i] = !string_i[i]

                # Convert the string into real values.
                vars = bitstrtonum(mgeod.design_vars, string_i, factors)

                # Evaluate the objective functions.
                (valid, f) = f_obj(vars)

                # Check if the solution is valid.
                if valid
                    # Create the candidate point.
                    candidate_point =
                        Pareto_Point{Nv, Nf}(vars, SVector{Nf, Float64}(f))

                    # Add the result to the rank.
                    @inbounds f_rank[i] = sRank(true, i, f[chosen_func])
                    @inbounds candidate_points[i] = candidate_point
                else
                    @inbounds f_rank[i] = sRank(false, i, 0.0)
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
            #
            # In this case (MGEO Var), we fill flip 1 bit per design variable.

            @inbounds for i = 1:Nv
                # Get the position of the i-th design variable in the string.
                bit_i = mgeod.design_vars[i].index
                bit_f = bit_i + mgeod.design_vars[i].bits - 1

                # Get the rank relative to the bit flip of the i-th design
                # variable.
                f_rank_dv = @view f_rank[bit_i:bit_f]

                # Sort the rank.
                sort!(f_rank_dv, lt=(a,b)->begin
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

                # Choose a bit to be changed for the new generation.
                bit_accepted = false

                while !bit_accepted
                    # Sample a valid bit in the string considering the current
                    # design variable.
                    b_sample = rand(bit_i:bit_f) - bit_i + 1

                    # Accept the change with probability r^(-τ), where r is the
                    # rank of the bit.
                    Pk = b_sample^(-mgeod.τ)

                    if rand() <= Pk
                        # The bit is accepted, then exit the loop.
                        string[f_rank_dv[b_sample].index] =
                            !string[f_rank_dv[b_sample].index]
                        bit_accepted = true
                    end
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

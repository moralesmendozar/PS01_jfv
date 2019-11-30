module exercise7

export multi_grid

# ---------
# Packages
# ---------
using Parameters, LinearAlgebra, Interpolations, NLsolve

# ---------
#  Modules
# ---------
include("./src/gridmake.jl")
using .exercise7aux

function multi_grid(economy, steadyState, nMidPoints::Array, tValueFunctionInitial::Array)

    # ---------------
    # 1. Housekeeping
    # ---------------
    @unpack vGridZ, mTransitionZ, vGridA, mTransitionA = economy
    @unpack capitalSteadyState, laborOneSteadyState, laborTwoSteadyState= steadyState

    nGridZ = length(vGridZ)
    nGridA = length(vGridA)

    dicLabor = Dict()
    dicValueFunction = Dict()

    # ----------------------------
    # 2. Auxiliary functions
    # ---------------------------
    # 2.1.
    function labor_choice(gridCapital, gridCapitalNext, gridProductivityZ, gridProductivityA, guess1, guess2, economy)

        @unpack α, β, δ = economy
        eZ = exp(gridProductivityZ)

        function focs12!(F,x)
             c1 = eZ * gridCapital^α * x[1]^(1-α) + (1-δ)*gridCapital - gridCapitalNext

             # FOC w.r.t labor 1
             F[1] = θ * (1-α)*eZ*gridCapital^α*x[1]^(-α) * (c1)^(θ-1)*(gridProductivityA*x[2])^(1-θ) - x[1] - x[2]

             # FOC w.r.t. labor 2
             F[2] = ( 1 - θ)* (c1)^θ *(gridProductivityA*x[2])^(-θ)*gridProductivityA - x[1] - x[2]
         end

        results = nlsolve(foc12!, [guess1, guess2])
        laborOne = results.zero[1]
        laborTwo = results.zero[2]

        return laborOne, laborTwo
    end

    # 2.2.
    function value_provisional(gridCapital::Real, gridCapitalNext::Real, gridProductivityZ::Real, gridProductivityA::Real, laborOne::Real, laborTwo::Real, expected::Real, economy)

        # 1. Housekeeping
        @unpack α, β, δ, θ = economy

        # 2. Consumption
        consumptionOne = exp(gridProductivityZ)*gridCapital^α * laborOne^(1-α) + (1-δ)*gridCapital - gridCapitalNext

        consumptionTwo = gridProductivityA * laborTwo

        # 3. Value provisional
        if consumptionOne <= 0 || consumptionTwo <= 0
            valueProvisional = -10000.0
        else
            utility = consumptionOne^(θ)*consumptionTwo^(1-θ) - 0.5*(laborOne+laborTwo)^2

            valueProvisional = (1-β)*utility + β*expected
        end

        return valueProvisional
    end

    # --------------------
    # 3. Multigrid scheme
    # --------------------

    for iGrid in 1:length(nMidPoints)

        # ----------------
        # 3.1. Grid Capital
        # ----------------
        vGridCapital = exericise7aux.gridmake(capitalMinimum, capitalMaximum, nGridCapital[iGrid])
        nGridCapital = length(vGridCapital)

        # -----------------
        # 3.2. Pre-allocation
        # -----------------
        tValueFunction = tValueFunctionInitial = fill(villaverdeSteadyState.utilitySteadyState, nGridCapital, nGridZ, nGridA)
        tValueFunctionNew = zeros(nGridCapital, nGridZ, nGridA)
        tPolicyFunctionIndex = zeros(Int64,nGridCapital,nGridZ, nGridA)
        tLaborOnePolicy = zeros(nGridCapital, nGridZ, nGridA)
        tLaborTwoPolicy = zeros(nGridCapital, nGridZ, nGridA)
        d4ProvisionalLaborOne = zeros(nGridCapital, nGridCapital, nGridZ, nGridA)
        d4ProvisionalLaborTwo = zeros(nGridCapital, nGridCapital, nGridZ, nGridA)
        vProvisionalValues = zeros(nGridCapital)
        vMaxDifference = Float64[]

        # ------------------
        # 3.3. Labor choices
        # ------------------
        guess1 = laborOneSteadyState
        guess2 = laborTwoSteadyState

        println("labor choices for all k' ... ")
        println(" ")
        for iA in 1:nGridA

            for iZ in 1:nGridZ

                for iCapital in 1:nGridCapital

                    for iCapitalNext in 1:nGridCapital

                        # Optimal labor choices
                        keyLabor = (vGridCapital[iCapital], vGridCapital[iCapitalNext], vGridZ[iZ],vGridA[iA])

                        if keyLabor ∈ keys(dicLabor)

                            laborOne, laborTwo = dicLabor[key]

                        else

                            laborOne, laborTwo = try
                                labor_choice(vGridCapital[iCapital], vGridCapital[iCapitalNext], vGridZ[iZ],vGridA[iA], guess1, guess2)
                            catch
                                guess1, guess2
                            end

                            # store in dictionary
                            keyLabor = (vGridCapital[iCapital], vGridCapital[iCapitalNext], vGridZ[iZ],vGridA[iA])

                            dicLabor[keyLabor] = [laborOne, laborTwo]
                        end
                    end
                end
            end
        end

        # ----------------------------
        # 3.4. Value Function Iteration
        # ----------------------------
        maxDifference = 1.0
        tolerance = 1.0e-6
        iteration = 0


        println("VFI ...")
        println(" ")
        while(maxDifference > tolerance)

            for iA in 1:nGridA

                for iZ in 1:nGridZ

                    for iCapital in 1:nGridCapital

                        for iCapitalNext in 1:nGridCapital

                            expected = 0.0

                            for iAnext in 1:nGridA

                                for iZnext in 1:nGridZ

                                    expected = expected + mTransitionZ[iZ, iZnext] * mTransitionA[iA, iAnext] * tValueFunction[iCapitalNext,iZnext,iAnext]
                                end
                            end

                            # labor
                            keyLabor = (vGridCapital[iCapital], vGridCapital[iCapitalNext], vGridZ[iZ],vGridA[iA])

                            laborOne, laborTwo = dicLabor[keyLabor]

                            # Values for all k'
                            vProvisionalValues[iCapitalNext] = value_provisional(vGridCapital[iCapital], vGridCapital[iCapitalNext], vGridZ[iZ],vGridA[iA], laborOne, laborTwo, expected, economy)

                        end

                        keyValueFunction = (vGridCapital[iCapital], vGridZ[iZ], vGridA[iA])

                        if keyValueFunction ∈ keys(dicValueFunction)

                            tValueFunctionNew[iCapital, iZ, iA], tPolicyFunctionIndex[iCapital, iZ, iA] = dicValueFunction[keyValueFunction]

                        else

                            tValueFunctionNew[iCapital, iZ, iA], tPolicyFunctionIndex[iCapital, iZ, iA] = findmax(vProvisionalValues)

                            dicValueFunction[keyValueFunction] = [tValueFunctionNew[iCapital, iZ, iA], tPolicyFunctionIndex[iCapital, iZ, iA]]

                        end

                    end
                end
            end

            maxDifference = norm(tValueFunctionNew-tValueFunction)
            vMaxDifference = push!(vMaxDifference,maxDifference)

            tValueFunction, tValueFunctionNew = tValueFunctionNew, tValueFunction

            iteration = iteration+1
            if(mod(iteration,10)==0 || iteration == 1)
                println(" Iteration = ", iteration, " Sup Diff = ", maxDifference)
            end
        end

    end

    tPolicyFunction = vGridCapital[tPolicyFunctionIndex]

    return tValueFunction, tPolicyFunction, vGridCapital, vMaxDifference

end

end

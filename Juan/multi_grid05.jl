module exercise7

export multigrid

using Parameters, LinearAlgebra, Interpolations, NLsolve

function multigrid(economy, steadyState, scheme)

    # ---------------
    # 0. Housekeeping
    # ---------------
    @unpack vGridZ, mTransitionZ, vGridA, mTransitionA = economy
    @unpack capitalSteadyState, laborOneSteadyState, laborTwoSteadyState, utilitySteadyState = steadyState
    @unpack nGridCapitalFirst, nGridCapitalSecond, nGridCapitalLast, cutoffFirstGrid, cutoffSecondGrid = scheme

    nGridZ = length(vGridZ)
    nGridA = length(vGridA)

    # -----------------
    # 1. Pre-allocation
    # -----------------
    tValueFunction = fill(utilitySteadyState, nGridCapitalLast, nGridZ, nGridA)
    tValueFunctionNew = zeros(nGridCapitalLast, nGridZ, nGridA)
    tPolicyFunctionIndex = zeros(Int64,nGridCapitalLast,nGridZ, nGridA)
    d4ProvisionalLaborOne = zeros(nGridCapitalLast, nGridCapitalLast, nGridZ, nGridA)
    d4ProvisionalLaborTwo = zeros(nGridCapitalLast, nGridCapitalLast, nGridZ, nGridA)
    vProvisionalValues = zeros(nnGridCapitalLast)
    vMaxDifference = Float64[]

    # ----------------------------
    # 2. Auxiliary functions
    # ---------------------------
    # 2.1.
    function labor_choice(gridCapital, gridCapitalNext, gridProductivityZ, gridProductivityA, guess1, guess2, economy)

        @unpack α, β, δ = economy
        eZ = exp(gridProductivityZ)

        function focs12!(F,x)
             c1 = eZ * gridCapital^α * x[1]^(1-α) + (1-δ)*gridCapital - gridCapitalNext

             # FOC w.r.t labor
             F[1] = θ * (1-α)*eZ*gridCapital^α*x[1]^(-α) * (c1)^(θ-1)*(gridProductivityA*x[2])^(1-θ) - x[1] - x[2]

             # FOC [l2]
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

    # ------------------
    # 3. Labor choices
    # ------------------
    guess1 = laborOneSteadyState
    guess2 = laborTwoSteadyState

    vGridCapital = collect(range(0.7 * capitalSteadyState, 1.3 * capitalSteadyState, length = nGridCapitalLast))

    println("labor choices for all k' ... ")
    println(" ")
    for iA in 1:nGridA

        for iZ in 1:nGridZ

            for iCapital in 1:nGridCapital

                for iCapitalNext in 1:nGridCapital

                    # Optimal labor choices
                    laborOne, laborTwo = try
                        labor_choice(vGridCapital[iCapital], vGridCapital[iCapitalNext], vGridZ[iZ],vGridA[iA], guess1, guess2)
                    catch
                        guess1, guess2
                    end

                    d4ProvisionalLaborOne[iCapitalNext, iCapital, iZ, iA] = laborOne
                    d4ProvisionalLaborTwo[iCapitalNext, iCapital, iZ, iA] = laborTwo

                end
            end
        end
    end

    # ----------------------------
    # 4. Value Function Iteration
    # ----------------------------
    maxDifference = 10.0
    tolerance = 1.0e-6
    iteration = 0

    println("VFI ...")
    println(" ")
    while(maxDifference > tolerance)

        if iteration < cutoffFirstGrid

            vGridCapital = collect(range(0.7 * capitalSteadyState, 1.3 * capitalSteadyState, length = nGridCapitalFirst))

            nGridCapital = length(vGridCapital)

        elseif cutoffFirstGrid <= iteration < cuttoffSecondGrid

            vGridCapital = collect(range(0.7 * capitalSteadyState, 1.3 * capitalSteadyState, length = nGridCapitalSecond))

            nGridCapital = length(vGridCapital)
        else

            vGridCapital = collect(range(0.7 * capitalSteadyState, 1.3 * capitalSteadyState, length = nGridCapitalLast))

            nGridCapital = length(vGridCapital)
        end

        itp = interpolate((vGridCapital, vGridZ, vGridA), tValueFunction, Gridded(Linear()))
        tValueFunctionInterpolate = itp(vGridCapital, vGridZ, vGridA)

        for iA in 1:nGridA

            for iZ in 1:nGridZ

                for iCapital in 1:nGridCapital

                    for iCapitalNext in 1:nGridCapital

                        expected = 0.0

                        for iAnext in 1:nGridA

                            for iZnext in 1:nGridZ

                                expected = expected + mTransitionZ[iZ, iZnext] * mTransitionA[iA, iAnext] * tValueFunctionInterpolate[iCapitalNext,iZnext,iAnext]
                            end
                        end

                        laborOne = d4ProvisionalLaborOne[iCapitalNext, iCapital, iZ, iA]

                        laborTwo = d4ProvisionalLaborTwo[iCapitalNext, iCapital, iZ, iA]

                        # Values for all k'
                        vProvisionalValues[iCapitalNext] = value_provisional(vGridCapital[iCapital], vGridCapital[iCapitalNext], vGridZ[iZ],vGridA[iA], laborOne, laborTwo, expected, economy)

                    end

                    tValueFunctionNew[iCapital, iZ, iA], tPolicyFunctionIndex[iCapital, iZ, iA] = findmax(vProvisionalValues)

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

    tPolicyFunction = vGridCapital[tPolicyFunctionIndex]

    return tValueFunction, tPolicyFunction, vGridCapital, vMaxDifference
end

end

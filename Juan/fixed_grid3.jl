module exercise3

export fixed_grid

using Parameters, LinearAlgebra, Interpolations, NLsolve

function fixed_grid(economy, steadyState)

    # ---------------
    # 0. Housekeeping
    # ---------------
    @unpack vGridZ, mTransitionZ, vGridA, mTransitionA = economy
    @unpack capitalSteadyState, laborOneSteadyState, laborTwoSteadyState = steadyState

    nGridCapital = 250
    nGridZ = length(vGridZ)
    nGridA = length(vGridA)

    # initial guess for labor choices
    guess1 = laborOneSteadyState
    guess2 = laborTwoSteadyState

    # ----------------
    # 1. Grid Capital
    # ----------------
    vGridCapital = collect(range(0.7 * capitalSteadyState, 1.3 * capitalSteadyState, length = nGridCapital))

    # -----------------
    # 2. Pre-allocation
    # -----------------
    tValueFunction = zeros(nGridCapital,nGridZ, nGridA)
    tValueFunctionNew = zeros(nGridCapital,nGridZ, nGridA)
    tPolicyFunctionIndex = zeros(Int64,nGridCapital,nGridZ, nGridA)
    tLaborOnePolicy = fill(laborOneSteadyState, nGridCapital,nGridZ, nGridA)
    tLaborTwoPolicy = fill(laborTwoSteadyState, nGridCapital,nGridZ, nGridA)
    vProvisionalValues = zeros(nGridCapital)

    # ----------------------------
    # 3. Value function iteration
    # ---------------------------

    # Auxiliary function
    function labor_choice(gridCapital, gridCapitalNext, gridProductivityZ, gridProductivityA, guess1, guess2)

        function parameters!(F,x, gridCapital, gridCapitalNext, gridProductivityZ, gridProductivityA)
            # fixed parameters
            α = 0.33; β = 0.96; δ = 0.9

            # changing parameters
            eZ = exp(gridProductivityZ)

            c1 = eZ * gridCapital^α * x[1]^(1-α) + δ*gridCapital - gridCapitalNext

        end

        function f2!(F,x)
            # Housekeeping
            α = 0.33; β = 0.96; δ = 0.9
            eZ = exp(gridProductivityZ)
            c1 = eZ * gridCapital^α * x[1]^(1-α) + δ*gridCapital - gridCapitalNext

            # FOC w.r.t. labor 1
            F[1] = 0.5 * (1-α) * eZ * gridCapital^α * x[1]^(-α) * (c1)^(-0.5) * (gridProductivityA * x[2])^(0.5) - x[1] - x[2]

            # FOC w.r.t. labor 2
            F[2] = 0.5 * gridProductivityA * (c1)^(0.5) * (gridProductivityA * x[2])^(-0.5) - x[1] - x[2]
        end

        results = nlsolve(f2!, [guess1, guess2])
        laborOne = results.zero[1]
        laborTwo = results.zero[2]

        return laborOne, laborTwo
    end


    function value_provisional(gridCapital::Real, gridCapitalNext::Real, gridProductivityZ::Real, gridProductivityA::Real, guess1::Real, guess2::Real, expected::Real, economy)

        # 1. Housekeeping
        @unpack α, β, δ = economy

        # 2. Choices
        # labor
        laborOne, laborTwo = try
            labor_choice(gridCapital, gridCapitalNext, gridProductivityZ, gridProductivityA, guess1, guess2)
        catch
            guess1, guess2
        end

        # consumptions
        consumptionOne = exp(gridProductivityZ)*gridCapital^α * laborOne^(1-α) + δ * gridCapital - gridCapitalNext

        consumptionTwo = gridProductivityA * laborTwo

        # 3. Value
        if consumptionOne <= 0 || consumptionTwo <= 0
            valueProvisional = -10000.0
        else
            utility = consumptionOne^(0.5)*consumptionTwo^(0.5) - 0.5*(1+laborTwo)^2

            valueProvisional = (1-β)*utility + β*expected
        end

        return valueProvisional, laborOne, laborTwo
    end

    # Main iteration
    maxDifference = 10.0
    tolerance = 1.0e-6
    iteration = 0
    count = 0

    while(maxDifference > tolerance)

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

                        vProvisionalValues[iCapitalNext], laborTwo, laborOne = value_provisional(vGridCapital[iCapital], vGridCapital[iCapitalNext], vGridZ[iZ],vGridA[iA],guess1, guess2, expected, economy)

                        guess1 = laborOne
                        guess2 = laborTwo

                    end

                    tValueFunctionNew[iCapital, iZ, iA], tPolicyFunctionIndex[iCapital, iZ, iA] = findmax(vProvisionalValues)

                end
            end
        end

        maxDifference = norm(tValueFunctionNew-tValueFunction)
        tValueFunction, tValueFunctionNew = tValueFunctionNew, tValueFunction

        iteration = iteration+1
        if(mod(iteration,10)==0 || iteration == 1)
            println(" Iteration = ", iteration, " Sup Diff = ", maxDifference)
        end
    end

    tPolicyFunction = vGridCapital[tPolicyFunctionIndex]

    return tValueFunction, tPolicyFunction, vGridCapital
end


end

module ex3

export 03_03_fixed_grid

using Parameters, LinearAlgebra, Interpolations

function 03_03_fixed_grid(economy, capitalSteadyState::Real)

    # 0. Get (unpack) the parameters back
    @unpack vGridZ, vGridA, mTranstnZA = economy

    nK = 200
    nZ = length(vGridZ)
    nA = length(vGridA)
    # Grid Capital
    vGridCapital = collect(range(0.7 * capitalSteadyState, 1.3 * capitalSteadyState, length = nGridCapital))

    mValueFunction = zeros(nK,nZ*nA)
    mValueFunctionNew = zeros(nK,nZ*nA)
    mPolicyFunctionIndex = zeros(Int64,nK,nZ*nA)
    vProvisionalValues = zeros(nK)

    # ----------------------------
    # 3. Value function iteration
    # ---------------------------

    # Auxiliary function
    function value_provisional(gridCapital::Real, gridCapitalNext::Real, gridProductivityZ::Real, gridProductivityA::Real, expected::Real, economy)

        # 1. Housekeeping
        @unpack α, β, δ = economy
        eZ = exp(gridProductivityZ)

        # 2. Choices
        # labor
        laborTwo = (1/((1-α)*eZ)) * ( eZ + δ*gridCapital^(1-α) - gridCapitalNext * gridCapital^(-α))

        # consumptions
        consumptionOne = eZ*gridCapital + δ * gridCapital - gridCapitalNext

        consumptionTwo = gridProductivityA * laborTwo

        # 3. Value
        if consumptionOne < 0 || consumptionTwo < 0
            valueProvisional = -10000.0
        else
            utility = consumptionOne^(0.5)*consumptionTwo^(0.5) - 0.5*(1+laborTwo)^2

            valueProvisional = (1-β)*utility + β*expected
        end

        return valueProvisional
    end

    # Main iteration
    maxDifference = 10.0
    tolerance = 1.0e-6
    iteration = 0

    while(maxDifference > tolerance)

        itp = interpolate((vGridCapital, vGridZ, vGridA), tValueFunction, Gridded(Linear()))
        tValueFunctionInterpolate = itp(vGridCapital, vGridZ,vGridA)

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

                        vProvisionalValues[iCapitalNext] = value_provisional(vGridCapital[iCapital], vGridCapital[iCapitalNext], vGridZ[iZ],vGridA[iA], expected, economy)

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

module exercise4vfi

export egm_vfi_standard

# ---------
# Packages
# ---------
using Parameters, NLsolve

function egm_vfi_standard(tValueFunctionTilde::Array, tMarketResourcesEndogenous::Array, economy, steadyState, vGridCapital::Array)

    # ---------------
    # 0. Housekeeping
    # ---------------
    @unpack α, θ, δ, β, vGridZ, mTransitionZ, vGridA, mTransitionA = economy
    @unpack capitalSteadyState, laborOneSteadyState, laborTwoSteadyState, utilitySteadyState = steadyState

    nGridZ = length(vGridZ)
    nGridA = length(vGridA)
    nGridCapital = length(vGridCapital)


    # -------------------
    # 2. Pre-allocation
    # -------------------
    tValueFunctionTildeNew = zeros(nGridCapital, nGridZ, nGridA)
    tPolicyFunction = zeros(nGridCapital, nGridZ, nGridA)
    tLaborOnePolicy = zeros(nGridCapital, nGridZ, nGridA)
    tLaborTwoPolicy = zeros(nGridCapital, nGridZ, nGridA)
    tConsumptionOnePolicy = zeros(nGridCapital, nGridZ, nGridA)
    tConsumptionTwoPolicy = zeros(nGridCapital, nGridZ, nGridA)


    # -----------------------
    # 3. Auxiliary functions
    # -----------------------
    # 3.1.
    function labor_choice(capital::Real, capitalNext::Real, eZ::Real, gridA::Real, guess1::Real, guess2::Real, economy)

        @unpack α, θ, δ, vGridZ = economy

        function f_mr!(F,x)

            # foc w.r.t. labor 1
            c11 = (θ*(1-α)*capital^α *x[1]^(-α) * x[2]^(1-θ)/ (x[1] + x[2]))^(1/(1-θ))

            # foc w.r.t. labor 2
            c12 = ((x[1] + x[2]) / (gridA * (1-θ)))^(1/θ) * gridA * x[2]

            # market resources
            F[1] = c11 + capitalNext - eZ*capital^α * x[1]^(1-α) - (1-δ)*capital

            F[2] = c12 + capitalNext - eZ*capital^α * x[1]^(1-α) - (1-δ)*capital

        end

        function j_mr!(J,x)

            aux1 = (θ*(1-α)*x*capital^α *x[1]^(-α) * x[2]^(1-θ)/ (x[1] + x[2]))^(1/(1-θ))

            aux2 = θ*(1-α)*eZ*capital^α * ((-α * x[1]^(-α-1)) / (x[1]+x[2]) - x[1]^(-α))

            aux3 =  aux2 * (aux1 * aux1^(-1)) * (gridA * x[2])

            J[1,1] = (1/1-θ) * aux3 - eZ * capital^α * (1-α) * x[1]^(-α)

            J[1,2] = aux1 * gridA

            aux4 = ( (x[1]+x[2]) / (gridA*(1-θ)) )^(1/θ)

            aux5 = (gridA*x[2]) / (gridA * θ*(1-θ))

            J[2,1] = aux5 * (aux4 * aux4^(-1))- eZ * capital^α * (1-α) * x[1]^(-α)

            J[2,2] = aux5 * (aux4 * aux4^(-1)) + gridA * aux4

        end

        results = nlsolve(f_mr!, j_mr!, [guess1, guess2])
        laborOne = results.zero[1]
        laborTwo = results.zero[2]

        return laborOne, laborTwo
    end

    # --------------------
    #  4. VFI in one step
    # --------------------
    # Initial guess for labor choices
    guess1 = laborOneSteadyState
    guess2 = laborTwoSteadyState

    println("Standard VFI in 1 step ...")
    println(" ")
    for iA in 1:nGridA

        for iZ in 1:nGridZ

            eZ = exp(vGridZ[iZ])
            gridCapitalNextPeriod = 1

            for iCapital in 1:nGridCapital

                valueHighSoFar = -10000.0
                capitalChoice  = vGridCapital[1]

                for iCapitalNext in gridCapitalNextPeriod:nGridCapital

                    # optimal labors
                    laborOne, laborTwo = try
                        labor_choice(vGridCapital[iCapital], vGridCapital[iCapitalNext], eZ, vGridA[iA], guess1,  guess2, economy)
                    catch
                        guess1, guess2
                    end


                    # consumptions
                    consumptionOne = exp(vGridZ[iZ])*vGridCapital[iCapital]^α * laborOne^(1-α) + (1-δ)*vGridCapital[iCapital]- vGridCapital[iCapitalNext]

                    consumptionTwo = vGridA[iA] * laborTwo

                    if consumptionOne < 0 || consumptionTwo < 0

                        valueProvisional = valueHighSoFar

                    else

                        expected = 0.0
                        for iAnext in 1:nGridA

                            for iZnext in 1:nGridZ

                                expected = expected + β * mTransitionZ[iZ,iZnext] * mTransitionA[iA, iAnext] * tValueFunctionTilde[iCapitalNext, iZnext, iAnext]

                            end
                        end

                        utility = consumptionOne^(θ)*consumptionTwo^(1-θ) - 0.5*(laborOne+laborTwo)^2

                        valueProvisional = (1-β) * utility + β * expected

                    end

                    if(valueProvisional>valueHighSoFar)

                        valueHighSoFar = valueProvisional
                        capitalChoice = vGridCapital[iCapitalNext]
                        gridCapitalNextPeriod = iCapitalNext

                        tLaborOnePolicy[iCapital, iZ, iA] = laborOne
                        tLaborTwoPolicy[iCapital, iZ, iA] = laborTwo

                    else

                        break # we have achieved the maximum

                    end

                end

                tValueFunctionTildeNew[iCapital, iZ, iA] = valueHighSoFar
                tPolicyFunction[iCapital, iZ, iA] = capitalChoice

            end
        end
    end


    maxDifference = maximum(abs.(tValueFunctionTildeNew - tValueFunctionTilde))

    tValueFunctionTilde, tValueFunctionTildeNew = tValueFunctionTildeNew, tValueFunctionTilde


    if maxDifference > 1e-6
        println("It should converge in 1 iteration! Something might be wrong. Check")
        @show maxDifference
    end


    # Consumption policy
    for iA in 1:nGridA
        for iZ in 1:nGridZ
            for iCapital in 1:nGridCapital

                tConsumptionOnePolicy[iCapital, iZ, iA] = exp(vGridZ[iZ])*vGridCapital[iCapital]^α * tLaborOnePolicy[iCapital, iZ, iA]^(1-α) + (1-δ)*vGridCapital[iCapital] - tPolicyFunction[iCapital, iZ, iA]

                tConsumptionTwoPolicy[iCapital, iZ, iA] = vGridA[iA]*tLaborTwoPolicy[iCapital, iZ, iA]

            end
        end
    end

return tPolicyFunction, tConsumptionOnePolicy, tConsumptionTwoPolicy, tLaborOnePolicy, tLaborTwoPolicy

end

end

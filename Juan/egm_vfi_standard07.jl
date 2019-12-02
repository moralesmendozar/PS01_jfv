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
    @unpack α, θ, δ, vGridZ, mTransitionZ, vGridA, mTransitionA = economy
    @unpack capitalSteadyState, laborOneSteadyState, laborTwoSteadyState, utilitySteadyState = steadyState

    nGridZ = length(vGridZ)
    nGridA = length(vGridA)
    nGridCapital = length(vGridCapital)

    # -----------------------
    # 1. Auxiliary functions
    # -----------------------

    function grid_capital(capitalNext::Real, eZ::Real, laborOne::Real, laborTwo::Real, economy)

        @unapck α, θ, δ, vGridZ = economy

        function f_mr!(F,x)

            c1 = (θ*(1-α)*x[1]^α *laborOne^(-α) * laborTwo^(1-θ)/ (laborOne + laborTwo))^(1/(1-θ))

            F[1] = c1 + capitalNext - eZ*x[1]^α - (1-δ)*x[1]

        end

        function j_mr!(J,x)
            aux1 = (θ*(1-α)*x[1]^α *laborOne^(-α) * laborTwo^(1-θ)/ (laborOne + laborTwo))^(θ/(1-θ))

            aux2 = θ * (1-α) * α * x[1]^(α - 1) * laborOne^(-α) * laborTwo^(1-θ)

            aux3 =  aux2/(laborOne + laborTwo)*(1-θ)

            J[1] = aux3 * aux1 - eZ * α * x[1]^(α-1) - (1-δ)
        end

        results = nlsolve(f_mr!, j_mr!, 0.5)
        capitalEndogenous = results.zero[1]

        return capitalEndogenous
    end

    # --------------------
    #  2. VFI in one step
    # --------------------

    for iA in 1:nGridA

        for iZ in 1:nGridZ

            eZ = exp(vGridZ[iZ])
            gridCapitalNextPeriod = 1

            for iCapital in 1:nGridCapital

                valueHighSoFar = -1000.0
                capitalChoice  = vGridCapital[1]

                for iCapitalNext in gridCapitalNextPeriod:nGridCapital

                    # endogenous capital grid
                    capitalEndgonous = gridcapital(vGridCapital[iCapitalNext], eZ, laborOneSteadyState, laborTwoSteadyState, economy)

                    # consumptions
                    consumptionOne = exp(vGridZ[iZ])*capitalEndogenous^α * laborOneSteadyState^(1-α) + (1-δ)*capitalEndogenous - vGridCapital[iCapitalNext]

                    consumptionTwo = vGridA[iA] * laborTwoSteadyState

                    if consumptionOne < 0
                        valueProvisional = valueHighSoFar
                    else
                        utility = consumptionOne^(θ)*consumptionTwo^(1-θ) - 0.5*(laborOneSteadyState+laborTwoSteadyState)^2

                        expected = 0.0
                        for iAnext in 1:nGridA

                            for iZnext in 1:nGridZ

                                expected = expected + β * mTransitionZ[iZ,iZnext] * mTransitionA[iA, iAnext] * tValueFunctionTilde[iCapitalNext, iZnext, iAnext]

                            end
                        end

                        valueProvisional = (1-β) * utility + β * expected
                    end

                    if(valueProvisional>valueHighSoFar)
                        valueHighSoFar = valueProvisional
                        capitalChoice = vGridCapital[iCapitalNext]
                        gridCapitalNextPeriod = iCapitalNext
                    else
                        break # We break when we have achieved the max
                    end

                end

                tValueFunctionNew[iCapital, iZ, iA] = valueHighSoFar
                tPolicyFunction[iCapital, iZ, iA] = capitalChoice
            end
        end
    end

    diff = norm(tValueFunctionTildeNew - tValueFunctionTilde)
    tValueFunctionTilde, tValueFunctionTildeNew = tValueFunctionTildeNew, tryValueFunctionTilde

    if diff > 1e-6
        error("It must converge in 1 iteration")
    end

    # Consumption policy
    for iA in 1:nGridA
        for iZ in 1:nGridZ
            for iCapital in 1:nGridCapital

                tConsumptionPolicy[iCapital, iZ, iA] = exp(vGridZ[iZ])*vGridCapital[iCapital]^α * laborOneSteadyState^(1-α) + (1-δ)*vGridCapital[iCapital] - tPolicyFunction[iCapitalNext, iA, iZ]

            end
        end
    end

return tPolicyFunction, tConsumptionPolicy

end

end

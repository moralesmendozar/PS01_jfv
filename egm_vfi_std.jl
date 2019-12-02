module ex4vfi

export egm_vfi_std

# Packages
using Parameters, NLsolve

function egm_vfi_std(tVFTilde::Array, tMarketResourcesEndogenous::Array, economy, steadyState, vGridK::Array)

    # Parameters:
    @unpack α,β, θ, δ, vGridZ, mTranstnZ, vGridA, mTranstnA = economy
    @unpack kss, l1ss, l2ss, utilitySS = steadyState

    nZ = length(vGridZ)
    nA = length(vGridA)
    nK = length(vGridK)

    # Auxiliary function with functions:

    function grid_capital(capitalNext::Real, eZ::Real, laborOne::Real, laborTwo::Real, economy)

        @unpack α, β, θ, δ, vGridZ = economy

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
        KEndog = results.zero[1]

        return KEndog
    end #function grid_capital(..)

    # VFI in one step
    tVFTildeNew = zeros(nK,nZ,nA)
    tPolicyFunction = zeros(nK,nZ,nA)
    for iA in 1:nA
        for iZ in 1:nZ

            eZ = exp(vGridZ[iZ])
            gridKNextPeriod = 1

            for iK in 1:nK

                valueHighSoFar = -1000.0
                KChoice  = vGridK[1]

                for iKNext in gridKNextPeriod:nK

                    # endogenous capital grid
                    KEndog = try
                        grid_capital(vGridK[iKNext], eZ, l1ss, l2ss, economy)
                    catch
                        vGridK[iKNext]
                    end

                    # consumptions
                    c1 = exp(vGridZ[iZ])*KEndog^α * l1ss^(1-α) + (1-δ)*KEndog - vGridK[iKNext]
                    c2 = vGridA[iA] * l2ss

                    if c1 < 0
                        valueProvisional = valueHighSoFar
                    else
                        utility = c1^(θ)*c2^(1-θ) - 0.5*(l1ss+l2ss)^2

                        expected = 0.0
                        for iAnext in 1:nA
                            for iZnext in 1:nZ
                                expected = expected + β * mTranstnZ[iZ,iZnext] * mTranstnA[iA, iAnext] * tVFTilde[iKNext, iZnext, iAnext]
                            end
                        end

                        valueProvisional = (1-β) * utility + β * expected
                    end

                    if(valueProvisional>valueHighSoFar)
                        valueHighSoFar = valueProvisional
                        KChoice = vGridK[iKNext]
                        gridKNextPeriod = iKNext
                    else
                        break # We break when we have achieved the max
                    end

                end

                tVFTildeNew[iK, iZ, iA] = valueHighSoFar
                tPolicyFunction[iK, iZ, iA] = KChoice
            end
        end
    end

    #diff = norm(tVFTildeNew - tVFTilde)
    diff = maximum(abs.(tVFTildeNew - tVFTilde))
    tVFTilde, tVFTildeNew = tVFTildeNew, tVFTilde

    if diff > 1e-6
        #error("Error! :( It must converge in 1 iteration")
        println("Error! :( It must converge in 1 iteration")
    end

    tConsPolicy = zeros(nK,nZ,nA)
    # Consumption policy
    for iA in 1:nA
        for iZ in 1:nZ
            for iK in 1:nK

                tConsPolicy[iK, iZ, iA] = exp(vGridZ[iZ])*vGridK[iK]^α * l1ss^(1-α) + (1-δ)*vGridK[iK] - tPolicyFunction[iK, iZ,iA]
                #tConsPolicy[iK, iZ, iA] = exp(vGridZ[iZ])*vGridK[iK]^α * l1ss^(1-α) + (1-δ)*vGridK[iK] - tPolicyFunction[iKNext, iA, iZ]

            end
        end
    end

    #tPolicyFunction = vGridK[tPolicyFunctionIndex]

    return tPolicyFunction, tConsPolicy

end

end

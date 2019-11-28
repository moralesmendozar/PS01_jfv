module ex3
export a_fixed_grid
using Parameters, LinearAlgebra, Interpolations

function a_fixed_grid(economy, steadyStateValues,nK)

    # 0. Get (unpack) the parameters back
    @unpack vGridZ, vGridA, mTranstnZA = economy
    @unpack kss, l1ss, l2ss = steadyStateValues
    nZ = length(vGridZ)
    nA = length(vGridA)

    # Grid Capital
    vGridK = collect(range(0.7 * kss, 1.3 * kss, length = nK))

    # Pre-allocation
    #   Create empty matrices where to record results (to return also) (:)
    mVF = zeros(nK,nZ*nA)
    mVFNew = zeros(nK,nZ*nA)
    mPolicyFnIndx = zeros(Int64,nK,nZ*nA)
    vProvisional = zeros(nK)

    # Initial guess for labour choices
    guessL1 = l1ss
    guessL2 = l2ss








    # Auxiliary function
    function val_provisional(gridK::Real, gridKNext::Real, gridProdZ::Real, gridProdA::Real, eV::Real, econ)
        # Get parameters:
        @unpack α, β, δ = econ
        eZ = exp(gridProdZ)

        # Labor:
        l2 = (1/((1-α)*eZ)) * ( eZ + δ*gridK^(1-α) - gridKNext * gridK^(-α))

        # consumptions
        c1 = eZ*gridK + δ * gridK - gridKNext

        c2 = gridProdA * l2

        # 3. Value
        if consumptionOne < 0 || consumptionTwo < 0
            valueProvisional = -10000.0
        else
            utility = consumptionOne^(0.5)*consumptionTwo^(0.5) - 0.5*(1+laborTwo)^2

            valueProvisional = (1-β)*utility + β*eV
        end

        return valueProvisional
    end






    # VFI
    maxDiff = 10.0
    tol = 1.0e-6
    iteration = 0

    while(maxDiff > tol)

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

                        vProvisionalValues[iCapitalNext] = val_provisional(vGridCapital[iCapital], vGridCapital[iCapitalNext], vGridZ[iZ],vGridA[iA], expected, economy)

                    end

                    mVFNew[iCapital, iZ, iA], mPolicyFnIndx[iCapital, iZ, iA] = findmax(vProvisionalValues)
                end
            end
        end

        maxDiff = norm(mVFNew-mVF)
        mVF, mVFNew = mVFNew, mVF

        iteration = iteration+1
        if(mod(iteration,10)==0 || iteration == 1)
            println(" Iteration = ", iteration, " Sup Diff = ", maxDiff)
        end
    end

    mPolicyFn = vGridK[mPolicyFnIndx]
    return mVF, mPolicyFn, vGridK

end

end

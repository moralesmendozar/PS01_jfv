module exercise4steadystate

export egm_labor_ss

# ---------
# Packages
# ---------
using Parameters
using LinearAlgebra, Interpolations, NLsolve

function egm_labor_ss(economy, steadyState, vGridCapitalNext::Array, vInitialGuess::Array)

    # ---------------
    # 0. Housekeeping
    # ---------------
    @unpack α, θ, δ, vGridZ, mTransitionZ, vGridA, mTransitionA = economy
    @unpack capitalSteadyState, laborOneSteadyState, laborTwoSteadyState, utilitySteadyState = steadyState

    nGridZ = length(vGridZ)
    nGridA = length(vGridA)
    nGridCapital = length(vGridCapitalNext)

    step = vGridCapitalNext[2]-vGridCapitalNext[1] # for derivatives

    # --------------------------
    # 2. Grid Market Resources
    # --------------------------
    eZ = exp.(vGridZ)
    productionOne = vGridCapitalNext.^(α) .* laborOneSteadyState^(1-α)

    mMarketResourcesNext = (ones(nGridCapital) * eZ') .* (productionOne * ones(1, nGridZ)) + (1-δ) .* (vGridCapitalNext * ones(1,nGridZ))

    tMarketResourcesNext = repeat(mMarketResourcesNext, 1, nGridA)

    # -----------------
    # 3. Pre-allocation
    # -----------------
    tValueFunctionTilde = repeat(vInitialGuess, 1, nGridZ, nGridA)
    tValueFunctionTildeNew = zeros(nGridCapital, nGridZ, nGridA)
    tDerivativeValueFunctionTilde = zeros(nGridCapital, nGridZ, nGridA)
    tOptimalConsumptionOne = zeros(nGridCapital, nGridZ, nGridA)
    tValueEndogenous = zeros(nGridCapital, nGridZ, nGridA)
    tMarketResourcesEndogenous = zeros(nGridCapital, nGridZ, nGridA)
    tGridCapitalEndogenous = zeros(nGridCapital, nGridZ, nGridA)

    # -----------------------------
    # 4. Value Function Iteration
    # -----------------------------
    maxDifference = 1.0
    tolerance = 1.0e-6
    iteration = 0

    println("VFI EGM fixed labor...")
    println(" ")
    while(maxDifference > tolerance)

        # 4.1. Derivatives at grid points only
        tDerivativeValueFunctionTilde[1,:,:] = (tValueFunctionTilde[2,:,:] - tValueFunctionTilde[1,:,:])/step

        tDerivativeValueFunctionTilde[end,:,:] = (tValueFunctionTilde[end,:,:] - tValueFunctionTilde[end-1,:,:])/step

        for iCapitalNext in 2:nGridCapital-1

            tDerivativeValueFunctionTilde[iCapitalNext,:,:] = (tValueFunctionTilde[iCapitalNext+1,:,:] - tValueFunctionTilde[iCapitalNext,:,:]) / step

        end

        # 4.2. Optimal consumption
        for iA in 1:nGridA

            tOptimalConsumptionOne[:,:,iA] = (tDerivativeValueFunctionTilde[:,:,iA]/θ).^(1/(θ-1)) * vGridA[iA] * laborTwoSteadyState

        end

        # 4.3. Update value function endogenous
        for iA in 1:nGridA

            tValueEndogenous[:,:,iA] = tOptimalConsumptionOne[:,:,iA].^θ * (vGridA[iA] * laborTwoSteadyState).^(1-θ) - fill( 0.5*(laborOneSteadyState + laborTwoSteadyState)^2, nGridCapital, nGridZ)+ tValueFunctionTilde[:,:,iA]

            for iZ in 1:nGridZ

                tMarketResourcesEndogenous[:,iZ, iA] = tOptimalConsumptionOne[:, iZ, iA] + vGridCapitalNext

            end
        end

        # 4.4. Interpolate on tomorrow's market resources grid ???
        int = interpolate((tMarketResourcesEndogenous[]), tValueEndogenous, Gridded(Linear()))

        tValueEndogenousNew = int((tMarketResourcesNext[]), tValueEndogenous)

        # 4.5. Compute the expectations
        for iA in 1:nGridA

            for iZ in 1:nGridZ

                for iAnext in 1:nGridA

                    for iZnext in 1:nGridZ

                        tValueFunctionTildeNew[:, iZ, iA] = tValueFunctionTildeNew[:, iZ, iA] + β * mTransitionZ[iZ,iZnext] * mTransitionA[iA, iAnext] * tValueEndogenousNew[:, iZnext, iAnext]

                    end
                end
            end
        end

        maxDifference = norm(tValueFunctionTildeNew-tValueFunctionTilde)

        tValueFunctionTilde, tValueFunctionTildeNew = tValueFunctionTildeNew, tValueFunctionTilde

        iteration = iteration+1
        if(mod(iteration,10)==0 || iteration == 1)
            println(" Iteration = ", iteration, " Sup Diff = ", maxDifference)
        end
    end

    return tValueEndogenous, tMarketResourcesEndogenous, vGridCapital

end

end

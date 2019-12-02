module ex4ss

export egm_lss

# Packages:
using Parameters
using LinearAlgebra, Interpolations, NLsolve

function egm_lss(economy, steadyState, vGridKNext::Array, vInitialGuess::Array)

    # Parameters
    @unpack α, β, θ, δ, vGridZ, mTranstnZ, vGridA, mTranstnA = economy
    @unpack kss, l1ss, l2ss, utilitySS = steadyState
    nZ = length(vGridZ)
    nA = length(vGridA)
    nK = length(vGridKNext)
    step = vGridKNext[2]-vGridKNext[1] # for derivatives

    # Grid Market Resources
    eZ = exp.(vGridZ)
    productionOne = vGridKNext.^(α) .* l1ss^(1-α)

    mMarketResourcesNext = (ones(nK) * eZ') .* (productionOne * ones(1, nZ)) + (1-δ) .* (vGridKNext * ones(1,nZ))
    tMarketResourcesNext = repeat(mMarketResourcesNext, 1, 1, nA)

    # PreAllocation
    tVFTilde = repeat(vInitialGuess, 1, nZ, nA)
    tVFTildeNew = zeros(nK, nZ, nA)
    tDeriVFTilde = zeros(nK, nZ, nA)
    tOptCons1 = zeros(nK, nZ, nA)
    tValueEndogenous = zeros(nK, nZ, nA)
    tMarketResourcesEndogenous = zeros(nK, nZ, nA)
    tGridKEndog = zeros(nK, nZ, nA)

    # VFI:
    maxDifference = 1.0
    tol = 1.0e-7
    iteration = 0

    println("VFI EGM fixed labor...")
    println(" ")

    while(maxDifference > tol)

        # 4.1. Derivatives at grid points only
        tDeriVFTilde[1,:,:] = (tVFTilde[2,:,:] - tVFTilde[1,:,:])/step
        tDeriVFTilde[end,:,:] = (tVFTilde[end,:,:] - tVFTilde[end-1,:,:])/step

        for iKNext in 2:nK-1
            tDeriVFTilde[iKNext,:,:] = (tVFTilde[iKNext+1,:,:] - tVFTilde[iKNext,:,:]) / step
        end

        # 4.2. Optimal consumption
        for iA in 1:nA
            tOptCons1[:,:,iA] = (tDeriVFTilde[:,:,iA]./θ).^(1/(θ-1)) * vGridA[iA] * l2ss
        end

        # 4.3. Update value function endogenous
        for iA in 1:nA
            #tValueEndogenous[:,:,iA] = tOptCons1[:,:,iA].^θ * (vGridA[iA] * l2ss).^(1-θ) - fill( 0.5*(l1ss + l2ss)^2, nK, nZ)+ tVFTilde[:,:,iA]
            tValueEndogenous[:,:,iA] = tOptCons1[:,:,iA].^θ * (vGridA[iA] * l2ss).^(1-θ) - fill( 0.5*(l1ss + l2ss)^2, nK, nZ)+ tVFTilde[:,:,iA]
            for iZ in 1:nZ
                tMarketResourcesEndogenous[:,iZ, iA] = tOptCons1[:, iZ, iA] + vGridKNext
            end
        end

        # 4.4. Interpolate on tomorrow's market resources grid ????
        tValueEndogenousNew = zeros(nK,nZ,nA)
        for iZ in 1:nZ
            for iA in 1:nA
                #itpi = interpolate(tMarketResourcesEndogenous[:,iZ,iA], tValueEndogenous[:,iZ,iA], Gridded(Linear()))
                itpi = LinearInterpolation(tMarketResourcesEndogenous[:,iZ,iA], tValueEndogenous[:,iZ,iA]; extrapolation_bc=Line())
                vectorTemp = itpi(tMarketResourcesNext[:,iZ,iA])
                tValueEndogenousNew[:,iZ,iA] = vectorTemp
            end
        end

        tVFTildeNew = zeros(nK, nZ, nA)
        # 4.5. Compute the expectations
        for iA in 1:nA
            for iZ in 1:nZ
                for iAnext in 1:nA
                    for iZnext in 1:nZ
                        tVFTildeNew[:, iZ, iA] = tVFTildeNew[:, iZ, iA] + β * mTranstnZ[iZ,iZnext] * mTranstnA[iA, iAnext] * tValueEndogenousNew[:, iZnext, iAnext]
                    end
                end
            end
        end

        maxDifference = norm(tVFTildeNew-tVFTilde)
        tVFTilde, tVFTildeNew = tVFTildeNew, tVFTilde

        iteration = iteration+1
        if(mod(iteration,10)==0 || iteration == 1)
            println(" Iteration = ", iteration, " Sup Diff = ", maxDifference)
        end
    end #while diff>tol

    #return tValueEndogenous, tMarketResourcesEndogenous
    return tVFTilde, tMarketResourcesEndogenous

end

end
